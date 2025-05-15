import os
import pandas as pd
import adaptive_sampling as msas

def scaling_surr():
    # Load data and pre-process by dropping irrelevant columns
    cwd = os.getcwd()
    filename = os.path.join(cwd, r'trial_area/Adaptive_sampling/humid_air_props_data.csv')
    ms_data = pd.read_csv(filename)
    ms_data.drop(['x_w_gas (mol/mol)', 'x_w_liq (mol/mol)', 'h_liq (J/mol)', 's_liq (J/mol/K)','v_liq (m3/mol)','x_gas (mol/mol)','RH','T_WB (K)',], axis=1, inplace=True)

    # Run adaptive sampling
    out_var = ['x_w_sat(mol/hmol)', 'h_gas (J/mol)', 's_gas (J/mol/K)','v_gas (m3/mol)']
    input_labels = ['T_DB (K)', 'P (Pa)', 'x_w (mol/mol)']
    input_bounds = {'T_DB (K)': [240, 400], 'P (Pa)': [50000, 200000], 'x_w (mol/mol)': [7.7e-7, 0.19]}

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
    filename = os.path.join(cwd, r'trial_area/Adaptive_sampling/humid_air_props_data.csv')
    HA_data = pd.read_csv(filename)
    HA_data.drop(['x_w_gas (mol/mol)', 'x_w_liq (mol/mol)', 'h_liq (J/mol)', 's_liq (J/mol/K)','v_liq (m3/mol)','x_gas (mol/mol)','RH','T_WB (K)',], axis=1, inplace=True)

    # Run adaptive sampling
    input_labels = ['T_DB (K)', 'P (Pa)', 'x_w (mol/mol)']
    input_bounds = {'T_DB (K)': [240, 400], 'P (Pa)': [50000, 200000], 'x_w (mol/mol)': [7.7e-7, 0.19]}

    # Initial sampling
    initial_no_samples = 50
    samples = msas.initial_sampling(HA_data, initial_no_samples, input_labels)

    # Model training
    out_vars = ['x_w_sat(mol/hmol)', 'h_gas (J/mol)', 's_gas (J/mol/K)','v_gas (m3/mol)']
    final_no_samples_max = 100
    model_type = 'gaussian'
    input_bounds = {'T_DB (K)': [240, 400], 'P (Pa)': [50000, 200000], 'x_w (mol/mol)': [7.7e-7, 0.19]}

    counter = 0
    while final_no_samples_max >= samples.shape[0]:
        print('\nIteration number: ', counter)
        print('=======================\n')
        model = msas.train_model_rbf(training_data=samples, basis_function=model_type, bounds=input_bounds, out_vars=out_vars)
        worst_surrogate = msas.compute_metrics(model=model, data=HA_data, out_var=out_vars)
        if final_no_samples_max != samples.shape[0]: # only update if required number of samples have not been selected.
            max_error_index = msas.determine_additional_point(model=model, data=HA_data, existing=samples, out_var=out_vars,worst_fit=worst_surrogate)
            samples = pd.concat([samples, HA_data.iloc[[max_error_index], :]])
            print('New selected point on iteration ', counter, ':', max_error_index)
            counter += 1
        else:
            break

    return model, samples


# Run mineral scaling example
model,samples = ph_example()
model.save_to_file('rbf_HA_adaptive_100.json', overwrite=True)
print(samples)

# # Run and plot pH example case to see newly added points:
# mod, samples = ph_example()