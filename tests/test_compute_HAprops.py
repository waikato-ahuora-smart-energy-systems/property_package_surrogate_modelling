# from compute_HAprops import (
#     compute_humid_air_properties,
#     generate_random_inputs,
#     generate_HAprops_dataset
# )

# import pandas as pd


# def test_valid_point():
#     result = compute_humid_air_properties(300, 101325, 0.01)
#     assert result is not None
#     expected_keys = [
#         "T_DB (K)", "P (Pa)", "x_w (mol/mol)", "x_w_gas (mol/mol)",
#         "h_gas (J/mol)", "s_gas (J/mol/K)", "v_gas (m3/mol)",
#         "x_w_liq (mol/mol)", "h_liq (J/mol)", "s_liq (J/mol/K)", "v_liq (m3/mol)",
#         "x_gas (mol/mol)", "RH", "T_WB (K)"
#     ]
#     for key in expected_keys:
#         assert key in result, f"Missing key: {key}"


# def test_invalid_T_above_sat():
#     result = compute_humid_air_properties(450, 101325, 0.01)
#     assert result is None


# def test_oversaturated_Y_handled():
#     result = compute_humid_air_properties(300, 101325, 0.5)
#     assert result is not None
#     assert result["RH"] == 1.0
#     assert result["T_WB (K)"] == 300
#     assert result["x_gas (mol/mol)"] <= 1.0


# def test_high_RH_edge_case():
#     result = compute_humid_air_properties(300, 101325, 0.035)
#     assert result is not None
#     assert 0.98 <= result["RH"] <= 1.0


# def test_valid_point_values_reasonable():
#     result = compute_humid_air_properties(300, 101325, 0.01)
#     assert result is not None
#     assert 250 < result["T_DB (K)"] < 350
#     assert 90000 < result["P (Pa)"] < 110000
#     assert 0 <= result["x_w (mol/mol)"] <= 0.5
#     assert 0 <= result["x_w_gas (mol/mol)"] <= 0.5
#     assert 0 < result["h_gas (J/mol)"] < 100000
#     assert 0 < result["s_gas (J/mol/K)"] < 500
#     assert 0 < result["v_gas (m3/mol)"] < 1
#     assert 0 <= result["RH"] <= 1
#     assert result["T_WB (K)"] < result["T_DB (K)"]


# def test_generate_random_inputs():
#     size = 10
#     T_bounds = (-70, 400)
#     P_bounds = (50, 200)
#     Y_bounds = (0.0, 0.5)

#     inputs = generate_random_inputs(size, T_bounds, P_bounds, Y_bounds)

#     assert len(inputs) == size
#     for T, P, Y in inputs:
#         assert 203.15 <= T <= 673.15
#         assert 50000 <= P <= 200000
#         assert 0.0 <= Y <= 0.5

# test_generate_random_inputs()

# def test_generate_HAprops_dataset_functional_serial():
#     df = generate_HAprops_dataset(
#         project_name="test_serial_run",
#         size=10,
#         T_bounds=(0, 50),
#         P_bounds=(90, 110),
#         Y_bounds=(0.01, 0.03),
#         save_csv=False,
#         serial=True
#     )

#     assert isinstance(df, pd.DataFrame)
#     assert not df.empty
#     expected_columns = [
#         "T_DB (K)", "P (Pa)", "x_w (mol/mol)", "x_w_gas (mol/mol)",
#         "h_gas (J/mol)", "s_gas (J/mol/K)", "v_gas (m3/mol)",
#         "x_w_liq (mol/mol)", "h_liq (J/mol)", "s_liq (J/mol/K)", "v_liq (m3/mol)",
#         "x_gas (mol/mol)", "RH", "T_WB (K)"
#     ]
#     for col in expected_columns:
#         assert col in df.columns
