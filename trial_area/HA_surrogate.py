import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# from keras import layers,regularizers
# from keras.layers import Layer
# from keras import backend as K
# import keras

from sklearn.metrics import mean_squared_error, r2_score
from sklearn import svm, tree
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.gaussian_process import kernels
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing, model_selection

from sklearn.metrics import r2_score
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.integrate import simpson

from idaes.core.surrogate.pysmo import polynomial_regression, radial_basis_function
from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D, surrogate_parity, surrogate_residual,
)


data_read = pd.read_csv(r'c:/Users/bjl25/Documents/property_package_surrogate_modelling/trial_area/humid_air_props_data.csv')

input_cols = ['T_DB (K)','P (Pa)','x_w (mol/mol)']

output_cols = ['x_w_sat(mol/hmol)','h_gas (J/mol)','s_gas (J/mol/K)','v_gas (m3/mol)']
# output_cols = ['Work_Mechanical']

inputs = data_read[input_cols]
outputs = data_read[output_cols]

# Randomly sample 100 data points from the dataset
sampled_data = data_read.sample(n=1000, random_state=42)

# Extract inputs and outputs from the sampled data
inputs = sampled_data[input_cols]
outputs = sampled_data[output_cols]

# Looking for points that are NaN in outputs
na_inx =  np.array(~outputs.isna().any(axis = 1))

inputs = inputs[na_inx] # must delete corresponding inputs as well!
outputs = outputs[na_inx]

# Train function takes the model object and input and output data for the training dataset
# Returns a trained model
def train_model(input_cols, output_cols, train_df, verbose = True):
  # Create the RBF trainer object
  rbf_trainer = PysmoRBFTrainer(input_labels=input_cols, output_labels=output_cols, training_dataframe = train_df)
  rbf_trainer.config.basis_function = 'gaussian'
  rbf_train = rbf_trainer.train_surrogate()
  surr = PysmoSurrogate(rbf_train, input_cols, output_cols)

  return surr

# Test function takes the model object and input data for the testing dataset and the names of the output variables so they are not mislabeled
# Returns a output model predictions for the testing dataset
def test_model(model, test_features):
  model_output = model.evaluate_surrogate(test_features)
  return model_output

def build_model():
  model = None
  return model


input_scaler = preprocessing.StandardScaler().fit(inputs)
output_scaler = preprocessing.StandardScaler().fit(outputs)

inputs_scaled = input_scaler.transform(inputs)
outputs_scaled = output_scaler.transform(outputs)

# Splitting into training and testing.
# For surrogate modleling it matters little but for real data it is crucial to think it through.

train_inputs, test_inputs = model_selection.train_test_split(inputs_scaled, test_size = 0.8, random_state=5)
train_outputs, test_outputs = model_selection.train_test_split(outputs_scaled, test_size = 0.8, random_state=5)
train_combined = pd.DataFrame(columns = input_cols+ output_cols, data = np.concatenate([train_inputs,train_outputs], axis = 1))
test_inputs_df = pd.DataFrame(columns = input_cols,data = test_inputs)


model = train_model(input_cols, output_cols, train_combined)
model_output = test_model(model, test_inputs_df) 


model_output = output_scaler.inverse_transform(model_output)
test_inputs = input_scaler.inverse_transform(test_inputs)
test_outputs = output_scaler.inverse_transform(test_outputs)
train_inputs = input_scaler.inverse_transform(train_inputs)
train_outputs = output_scaler.inverse_transform(train_outputs)

# for PYSMO

df_model_output = pd.DataFrame(model_output, columns = outputs.columns)
df_test_inputs = pd.DataFrame(test_inputs, columns = inputs.columns)
df_train_inputs = pd.DataFrame(train_inputs, columns = inputs.columns)
df_test_outputs = pd.DataFrame(test_outputs, columns = outputs.columns)
df_train_outputs = pd.DataFrame(train_outputs, columns = outputs.columns) 

full_data = np.concatenate([train_inputs,train_outputs],axis=1)
full_data.shape

av = np.average([np.average(abs(df_test_outputs[j] - df_model_output[j]))/np.average(df_model_output[j]) for j in df_model_output.columns])   
r2 = np.average([r2_score(df_test_outputs[j], df_model_output[j]) for j in df_model_output.columns]) 
av2 = (abs(df_test_outputs - df_model_output))/np.average(df_model_output)

print("Average Absolute Relative Error (av):", av)
print("Average R2 Score (r2):", r2)
print("Average Absolute Error per Output (av2):")
print(av2)


for i in output_cols:
    print(i)
    plt.figure()
    plt.title(i)
    plt.ylabel('model')
    plt.xlabel('data')
    plt.plot(df_test_outputs[i], df_model_output[i],'x')
    plt.show()

    