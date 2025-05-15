import numpy as np
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
import pandas as pd
import os

import matplotlib.pyplot as plt

# Load your surrogate model and true data

data_read = pd.read_csv(r'c:/Users/bjl25/Documents/property_package_surrogate_modelling/trial_area/humid_air_props_data_copy.csv')

input_cols = ['T_DB (K)','P (Pa)','x_w (mol/mol)']

output_cols = ['x_w_sat(mol/hmol)','h_gas (J/mol)','s_gas (J/mol/K)','v_gas (m3/mol)']

script_dir = os.path.dirname(__file__)
pysmo_surrogate = PysmoSurrogate.load_from_file(
            os.path.join(script_dir,"rbf_HA_10000.json")
        )

# Define test data 
test_data = data_read[input_cols]
# Define true outputs for the test data (
true_outputs = data_read[output_cols]

# Predict outputs using the surrogate model
predicted_outputs = pysmo_surrogate.evaluate_surrogate(test_data)

# Calculate absolute and percentage deviations
absolute_deviation = np.abs(predicted_outputs - true_outputs)
percentage_deviation = (absolute_deviation / true_outputs) * 100



# Plot absolute deviation for each output against temperature, pressure, and x_w
for i, col in enumerate(output_cols):
    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['T_DB (K)'], absolute_deviation[col], label=f'Absolute Deviation: {col} vs T_DB', alpha=0.7)
    plt.xlabel('Temperature (T_DB in K)')
    plt.ylabel('Absolute Deviation')
    plt.title(f'Absolute Deviation for {col} vs Temperature')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['P (Pa)'], absolute_deviation[col], label=f'Absolute Deviation: {col} vs P', alpha=0.7)
    plt.xlabel('Pressure (P in Pa)')
    plt.ylabel('Absolute Deviation')
    plt.title(f'Absolute Deviation for {col} vs Pressure')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['x_w (mol/mol)'], absolute_deviation[col], label=f'Absolute Deviation: {col} vs x_w', alpha=0.7)
    plt.xlabel('x_w (mol/mol)')
    plt.ylabel('Absolute Deviation')
    plt.title(f'Absolute Deviation for {col} vs x_w')
    plt.legend()
    plt.grid()
    plt.show()

# Plot percentage deviation for each output against temperature, pressure, and x_w
for i, col in enumerate(output_cols):
    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['T_DB (K)'], percentage_deviation[col], label=f'Percentage Deviation: {col} vs T_DB', alpha=0.7)
    plt.xlabel('Temperature (T_DB in K)')
    plt.ylabel('Percentage Deviation (%)')
    plt.title(f'Percentage Deviation for {col} vs Temperature')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['P (Pa)'], percentage_deviation[col], label=f'Percentage Deviation: {col} vs P', alpha=0.7)
    plt.xlabel('Pressure (P in Pa)')
    plt.ylabel('Percentage Deviation (%)')
    plt.title(f'Percentage Deviation for {col} vs Pressure')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.scatter(test_data['x_w (mol/mol)'], percentage_deviation[col], label=f'Percentage Deviation: {col} vs x_w', alpha=0.7)
    plt.xlabel('x_w (mol/mol)')
    plt.ylabel('Percentage Deviation (%)')
    plt.title(f'Percentage Deviation for {col} vs x_w')
    plt.legend()
    plt.grid()
    plt.show()

    # # Save the graphs
    # abs_dev_filename = os.path.join(script_dir, f'absolute_deviation_{col}.png')
    # plt.savefig(abs_dev_filename)
    # plt.close()

    # perc_dev_filename = os.path.join(script_dir, f'percentage_deviation_{col}.png')
    # plt.savefig(perc_dev_filename)
    # plt.close()