# import os
# import pytest
# import pandas as pd
# from train_surrogate import (
#     get_training_data_from_file,
#     train_surrogate_model,
#     get_model_validation_plots,
#     PysmoSurrogate
# )

# SAMPLE_ROWS = 50
# XMIN = [200, 50000, 0]
# XMAX = [600, 200000, 0.5]
# TEST_DATA_FILE = "tests/sample_humid_air_data.csv"
# N_INPUTS = 3
# TEMP_PATH = ".temp/"

# # ----------------------------------------------------------------------------------
# # Data loading

# def test_get_training_data():
#     df = get_training_data_from_file(TEST_DATA_FILE, sample_rows=SAMPLE_ROWS)
#     assert isinstance(df, pd.DataFrame)
#     assert df.shape[0] <= SAMPLE_ROWS
#     assert "temperature" in df.columns

# def test_file_not_found():
#     with pytest.raises(FileNotFoundError):
#         get_training_data_from_file("nonexistent.csv")

# # # ----------------------------------------------------------------------------------
# # # Training surrogate

# @pytest.mark.parametrize("model_type", ["poly", "rbf", "kriging"])
# def test_train_surrogate_model(model_type):
#     df = get_training_data_from_file(TEST_DATA_FILE, sample_rows=SAMPLE_ROWS)
    
#     model_path = f"{TEMP_PATH}test_model_{model_type}"

#     surrogate, train_data, val_data = train_surrogate_model(
#         df=df,
#         n_inputs=N_INPUTS,
#         xmin=XMIN,
#         xmax=XMAX,
#         model_type=model_type,
#         output_filename=str(model_path)
#     )

#     assert isinstance(surrogate, PysmoSurrogate)
#     assert not train_data.empty
#     assert not val_data.empty

# def test_invalid_model_type():
#     df = get_training_data_from_file(TEST_DATA_FILE, sample_rows=SAMPLE_ROWS)
#     with pytest.raises(ValueError):
#         train_surrogate_model(df, n_inputs=N_INPUTS, xmin=XMIN, xmax=XMAX, model_type="bad_model")

# def test_empty_dataframe():
#     with pytest.raises(ValueError):
#         train_surrogate_model(pd.DataFrame(), n_inputs=N_INPUTS, xmin=XMIN, xmax=XMAX, model_type="poly")

# # ----------------------------------------------------------------------------------
# # Plotting

# def test_model_validation_plots():
#     df = get_training_data_from_file(TEST_DATA_FILE, sample_rows=SAMPLE_ROWS)

#     model_path = f"{TEMP_PATH}plot_model"

#     try:
#         surrogate, train_data, _ = train_surrogate_model(
#             df,
#             n_inputs=N_INPUTS,
#             xmin=XMIN,
#             xmax=XMAX,
#             model_type="rbf",
#             output_filename=str(model_path)
#         )
#     except ValueError as e:
#         raise

#     get_model_validation_plots(surrogate, train_data, str(model_path), "train")
#     assert os.path.exists(f"results/{model_path}_train_scatter2D.pdf")
#     assert os.path.exists(f"results/{model_path}_train_parity.pdf")
#     assert os.path.exists(f"results/{model_path}_train_residual.pdf")

# def test_plotting_empty_data():
#     with pytest.raises(ValueError):
#         get_model_validation_plots(None, pd.DataFrame(), "dummy", "empty")
