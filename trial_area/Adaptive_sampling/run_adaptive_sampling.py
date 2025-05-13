import os
import pandas as pd
import adaptive_sampling as msas

def scaling_surr():
    # Load data and pre-process by dropping irrelevant columns
    cwd = os.getcwd()
    filename = os.path.join(cwd, r'')
    ms_data = pd.read_csv(filename)
    ms_data.drop(['ionicStrength', 'ionicStrengthXBased (mol/mol)', 'ionicStrengthMBased (mol/kg)', 'activityCoefficientsXBased', 'activityCoefficientsMBased'], axis=1, inplace=True)
    ms_data.rename(columns={'# Soda Ash Dose (mg/L)': 'NA2CO3',
                            'CO2 dose (mg/L)': 'CO2', 
                            'RO recovery': 'RO_recovery', 
                            'CASO4.2H2O': 'Gypsum'}, inplace=True)
    ms_data = ms_data[ms_data['RO_recovery'] <= 0.95]

    # Run adaptive sampling
    out_var = 'CACO3'
    input_labels = ['NA2CO3', 'RO_recovery', 'pH_pre', 'Pressure (atm)']
    input_bounds = {'NA2CO3': [0, 750], 'RO_recovery': [0.48, 0.95], 'pH_pre': [6, 9],  'Pressure (atm)': [10, 110]}

    model, ms_data_caco3_filter, ms_data_caco3_tightfilter = msas.scaling_adaptive_sampling_function(dataset=ms_data, 
                                                    out_var = out_var, 
                                                    input_labels = input_labels,
                                                    input_bounds = input_bounds, 
                                                    initial_no_samples=50, 
                                                    classification_adaptive_samples=50, 
                                                    mae_adaptive_samples=50, 
                                                    scaling_cutoff=2.0,
                                                    recovery_cutoff=0.95, 
                                                    tightfilter_low=0.75, 
                                                    tightfilter_high=1.25
                                                   )
    return model


def ph_example():
    # Load data
    cwd = os.getcwd()
    filename = os.path.join(cwd, r'data\brackish1_carbonation_from_scratch_extended.csv')
    carbonation_data = pd.read_csv(filename)

    # Initial sampling
    initial_no_samples = 50
    samples = msas.initial_sampling(carbonation_data, initial_no_samples, ['NA2CO3','CO2'])

    # Model training
    out_var = 'pH'
    final_no_samples_max = 100
    model_type = 'cubic'
    input_bounds = {'NA2CO3': [0, 1200], 'CO2': [0, 300]}

    counter = 0
    while final_no_samples_max >= samples.shape[0]:
        print('\nIteration number: ', counter)
        print('=======================\n')
        model = msas.train_model_rbf(training_data=samples, basis_function=model_type, bounds=input_bounds, out_var=out_var)
        msas.compute_metrics(model=model, data=carbonation_data, out_var=out_var)
        if final_no_samples_max != samples.shape[0]: # only update if required number of samples have not been selected.
            max_error_index = msas.determine_additional_point(model=model, data=carbonation_data, existing=samples, out_var=out_var)
            samples = pd.concat([samples, carbonation_data.iloc[[max_error_index], :]])
            print('New selected point on iteration ', counter, ':', max_error_index)
            counter += 1
        else:
            break

    return model, samples


# Run mineral scaling example
model = scaling_surr()

# # Run and plot pH example case to see newly added points:
# mod, samples = ph_example()