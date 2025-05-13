# Import statements
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Import IDAES libraries
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.alamopy import AlamoTrainer, AlamoSurrogate
from idaes.core.surrogate.pysmo import sampling as sp 
from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer, PysmoKrigingTrainer, PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.metrics import compute_fit_metrics
from idaes.core.surrogate.plotting.sm_plotter import surrogate_scatter2D, surrogate_parity, surrogate_residual
from idaes.core.surrogate.surrogate_block import SurrogateBlock


def initial_sampling(data, no_samples, x_labels):
	"""
	Selects initial samples with Hammersley sampling technique
	"""
	sp_load = sp.HammersleySampling(data, number_of_samples=no_samples, sampling_type='selection', xlabels=x_labels)
	samples = sp_load.sample_points()
	return samples


def train_model_rbf(training_data, basis_function, bounds, out_var=None): 
    # Create PySMO trainer object
    input_labels = list(bounds.keys()) 
    output_labels = ['CACO3'] if out_var is None else [out_var]
    carbonation_trainer_2 = PysmoRBFTrainer(input_labels=input_labels,
                            output_labels=output_labels,
                            training_dataframe = training_data)
    # Set PySMO options
    carbonation_trainer_2.config.basis_function = basis_function
    carbonation_trainer_2.config.regularization = True
  
    # Train surrogate (calls PySMO through IDAES Python wrapper)
    rbf_train = carbonation_trainer_2.train_surrogate()
    rbf_surr_carbonation = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds=bounds) 

    return rbf_surr_carbonation


def train_model_kriging(training_data, bounds, out_var=None): 
    # Create PySMO trainer object
    input_labels = list(bounds.keys()) 
    output_labels = ['CACO3'] if out_var is None else [out_var]
    carbonation_trainer_3 = PysmoKrigingTrainer(input_labels=input_labels,
                            output_labels=output_labels,
                            training_dataframe = training_data)

    # Set PySMO options
    carbonation_trainer_3.config.numerical_gradients = False
    carbonation_trainer_3.config.regularization = True
    

    # Train surrogate (calls PySMO through IDAES Python wrapper)
    krg_train = carbonation_trainer_3.train_surrogate()

    krg_surr_carbonation = PysmoSurrogate(krg_train, input_labels, output_labels, input_bounds=bounds) 

    return krg_surr_carbonation


def compute_metrics(model, data, out_var=None):
	"""
	Computes standard metrics
	"""
	var = 'CACO3' if out_var is None else out_var
	surrogate_models = [model]
	res_names = []
	for m in range(0, len(surrogate_models)):
	    err = compute_fit_metrics(surrogate_models[m], data)
	    err = pd.DataFrame.from_dict(err)
	    # print(err)
	    # res_names.append(err)
	    print('\nModel metrics: R2 =', err[var]['R2'], 'maxAE =', err[var]['maxAE'], '\n')


def determine_additional_point(model, data, existing, out_var=None):
	"""
	Determines and adds worst fit point (based on absolute error) to the training data.
	"""
	ms_data_mod = data.copy()
	# Drop columns already in training set - no point should exist in training set more than once.
	var = 'CACO3' if out_var is None else out_var
	unused_indexes = ms_data_mod.round(4).merge(existing.round(4),how='left',indicator=True).loc[lambda x : x['_merge']=='left_only'].index
	ms_data_mod = ms_data_mod.iloc[unused_indexes]
	pred_values = model.evaluate_surrogate(ms_data_mod)
	ms_data_mod['RBF'] = pred_values
	ms_data_mod['abs_error'] = abs(ms_data_mod['RBF'] - ms_data_mod[var])
	max_error_index = ms_data_mod['abs_error'].idxmax()
	print('Largest absolute deviation:', ms_data_mod['abs_error'].max())
	print('Row of interest:\n', ms_data_mod[ms_data_mod['abs_error']==ms_data_mod['abs_error'].max()])
	return max_error_index


def compute_ic_metrics(model, data, plot=False, out_var=None):
	var = 'CACO3' if out_var is None else out_var
	added_acid = 'CO2' if 'CO2' in data.columns else 'HCL' if 'HCL' in data.columns else 'H2SO4' if 'H2SO4' in data.columns else 'unknown'
	if added_acid == 'unknown':
		raise ValueError("Acid is not in expected list: ['HCL', 'CO2', 'H2SO4']")
	pred2 = model.evaluate_surrogate(data)
	xv = pd.DataFrame()
	xv['NA2CO3'] = data['NA2CO3'] 
	xv[added_acid] = data[added_acid] 
	xv['RO_recovery'] = data['RO_recovery']
	xv['pH_pre'] = data['pH_pre'] # xv['pH'] = data['pH']
	xv['Pressure (atm)'] = data['Pressure (atm)']
	xv['Actual ' + var] = data[var]
	xv[var] = data[var] # --- needed for other function
	xv['Predicted ' + var] = pred2[var]

	# Information criterion matrix
	print('Number of training points:', xv.shape[0])
	print('Number of true positives:', xv[(xv['Actual ' + var]> 1) & (xv['Predicted ' + var]> 1)].shape[0])
	print('Number of true negatives:', xv[(xv['Actual ' + var]< 1) & (xv['Predicted ' + var]< 1)].shape[0])
	print('Number of false positives:', xv[(xv['Actual ' + var]< 1) & (xv['Predicted ' + var]> 1)].shape[0])
	print('Number of false negatives:', xv[(xv['Actual ' + var]> 1) & (xv['Predicted ' + var]< 1)].shape[0])
	print('Fraction of correct predictions: ', (xv[(xv['Actual ' + var]> 1) & (xv['Predicted ' + var]> 1)].shape[0] +  xv[(xv['Actual ' + var]< 1) & (xv['Predicted ' + var]< 1)].shape[0]) / xv.shape[0])

	xv_correct = pd.concat(
	    [
	        xv[((xv['Actual ' + var]> 1) & (xv['Predicted ' + var]> 1))], 
	        xv[(xv['Actual ' + var]< 1) & (xv['Predicted ' + var]< 1) ]
	    ])

	xx = (pd.merge(xv, xv_correct, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1))

	if plot:
		xw = xv[(xv['Actual ' + var]> 0.5) & (xv['Actual ' + var]< 2)]
		plt.plot(xw['Actual ' + var], xw['Predicted ' + var], 'o')
		plt.plot(xw['Actual ' + var], xw['Actual ' + var])
		plt.xlim([0.4, 2.0])
		plt.show()

		print(xx)

	return xx


def scaling_adaptive_sampling_function(dataset, out_var, input_labels, input_bounds, initial_no_samples=50, classification_adaptive_samples=50, mae_adaptive_samples=50, scaling_cutoff=2.0, recovery_cutoff=1.0, tightfilter_low=0.0, tightfilter_high=100.0, initial_sampling_xcols=['NA2CO3','CO2','RO_recovery']):
    """
    Function for carrying out two-stage adaptive sampling for mineral scaling surrogate models.
    
    Samples are selected adaptively in two steps.
        - In the first step, ``classification_adaptive_samples`` points are selected to improve the classification accuracy. To do this, update points are only selected from the set of points that have been wrongly classified by the current model. The criterion for selecting from this filtered dataset is the maximum absolute error.
        - In the second step, ``mae_adaptive_samples`` points are selected to improve the MAE within a specified scaling tendency region. This sampling round allows us to improve the model performance within a specified, restricted ST range. The bounds of this region are defined by the ``tightfilter_low`` and ``tightfilter_high``. The criterion for selecting from this filtered dataset is the maximum absolute error.
    
    Args:
        dataset:				Input dataset containing independent variables and mineral scaling predictions 
        out_var: 				Output variable name in ``dataset``
        input_labels: 			List containing the variable names of the inputs in ``dataset``
        input_bounds: 			Dictionary containing the bounds of each input variable in ``input_labels``. 
        initial_number_samples: Initial number of samples selected via space-filling methods (Hammersley sampling)
        classification_adaptive_samples: Number of adaptive samples selected that focus on classification accuracy maximization. The points added to the training data here are selected from the list of points wrongly classified.
        mae_adaptive_samples: 	Number of samples selected purely based on MAE minimization only.
        scaling_cutoff: 	Threshold for scaling values in training dataset. Values higher than this are dropped from the dataset before model training starts. Default is 2.0.
        recovery_cutoff: 	Maximum recovery regime to be covered in training dataset. Data points with higher recovery values are dropped before model training starts. Default is 1.0.
        tightfilter_low: 	Lower scaling tendency bound for the restricted region dataset. Default is 0.
        tightfilter_high: 	Upper scaling tendency bound for the restricted region dataset. Default is 100.
        initial_sampling_xcols : List of column names for inputs to be used to select initial set of samples. Default is ['NA2CO3','CO2','RO_recovery'].
    
    Returns:
        model: Final model
        ms_data_filter: Full training dataset. Contains all potential points that were considered for selection after cutoff limits were applied.
        ms_data_tightfilter: Dataset of points within the restricted ST region specified by the user.
        
    """
    ms_data = dataset
    ms_data_filter = ms_data[ms_data[out_var] < scaling_cutoff]
    ms_data_filter = ms_data_filter[ms_data_filter['RO_recovery'] <= recovery_cutoff]
    ms_data_filter.reset_index(drop=True, inplace=True)

    # Create tight data with values supplied by user
    ms_data_tightfilter = ms_data_filter[(ms_data_filter[out_var] >= tightfilter_low) & (ms_data_filter[out_var] <= tightfilter_high)]

    # Initial sampling
    samples = initial_sampling(ms_data_filter, initial_no_samples, initial_sampling_xcols)

    # Model training
    final_no_samples_max = initial_no_samples + classification_adaptive_samples + mae_adaptive_samples
    model_type = 'cubic'

    counter = 0

    while final_no_samples_max >= samples.shape[0]:
        print('\nIteration number: ', counter)
        print('=======================\n')
        model = train_model_rbf(samples, model_type, input_bounds, out_var=out_var)
        compute_metrics(model, ms_data_filter, out_var=out_var)
        if final_no_samples_max != samples.shape[0]: # only update if required number of samples have not been selected.
            ms_data_failed = compute_ic_metrics(model, ms_data_filter, out_var=out_var)

            if samples.shape[0] < initial_no_samples + classification_adaptive_samples: # Check if all classification samples have been collected.
                try:
                    # Select additional point from dataframes of failed classification runs only. 
                    # Will fail if there are no failed classification runs or all failing runs are already in training data
                    max_error_index = determine_additional_point(model, ms_data_failed, samples, out_var=out_var)
                    print('Sample number', samples.shape[0], '; Classification selection.')
                except ValueError:
                    # If failure occurs due to above stated reasons, select from larger filtered dataset of < 3 
                    max_error_index = determine_additional_point(model, ms_data_filter, samples, out_var=out_var)
                    print('Sample number', samples.shape[0], '; Absolute error selection in classification group.')
            else:
                max_error_index = determine_additional_point(model, ms_data_tightfilter, samples, out_var=out_var)
                print('Sample number', samples.shape[0], '; Absolute error selection.')


            samples = pd.concat([samples, ms_data_filter.iloc[[max_error_index], :]])
            print('New selected point on iteration ', counter, ':', max_error_index)
            counter += 1
        else:
            break

    # Print final sample set
    model = train_model_rbf(samples, model_type, input_bounds, out_var=out_var)
    print(samples)
    compute_ic_metrics(model, ms_data, plot=True, out_var=out_var)
    compute_metrics(model, ms_data_filter, out_var=out_var)
    
    return model, ms_data_filter, ms_data_tightfilter