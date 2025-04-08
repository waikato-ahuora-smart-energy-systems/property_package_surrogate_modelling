import numpy as np
import pandas as pd
from CoolProp.CoolProp import HAPropsSI, PropsSI
from time import perf_counter as timer
from multiprocessing import Pool

# --- Constants ---
MW_A = PropsSI("M", "Air")
MW_W = PropsSI("M", "water")
Ttriple = PropsSI("Ttriple", "water")

# --- Core Property Function ---
def compute_humid_air_properties(T, P, Y):
    """
    Computes humid air thermodynamic properties at given dry-bulb temperature (T),
    pressure (P), and water mole fraction (Y).

    Returns a dictionary of properties if valid; otherwise, returns None.
    """
    try:
        T_w_sat_at_P = PropsSI("T", "P", P, "Q", 1, "water")
        if T >= T_w_sat_at_P:
            return None

        Y_sat = HAPropsSI("Y", "T", T, "P", P, "RH", 1)

        if Y <= Y_sat:
            MW_gas = Y * MW_W + (1 - Y) * MW_A
            RH = HAPropsSI("RH", "T", T, "P", P, "Y", Y)
            return {
                "T_DB (K)": T,
                "P (Pa)": P,
                "x_w (mol/mol)": Y,
                "x_w_gas (mol/mol)": Y,
                "h_gas (J/mol)": HAPropsSI("Hha", "T", T, "P", P, "Y", Y) * MW_gas,
                "s_gas (J/mol/K)": HAPropsSI("Sha", "T", T, "P", P, "Y", Y) * MW_gas,
                "v_gas (m3/mol)": HAPropsSI("Vha", "T", T, "P", P, "Y", Y) * MW_gas,
                "x_w_liq (mol/mol)": 1.0,
                "h_liq (J/mol)": PropsSI("H", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "s_liq (J/mol/K)": PropsSI("S", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "v_liq (m3/mol)": 1 / PropsSI("D", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "x_gas (mol/mol)": 1.0,
                "RH": RH,
                "T_WB (K)": HAPropsSI("Twb", "T", T, "P", P, "Y", Y)
            }
        else:
            x_w_gas = HAPropsSI("Y", "T", T, "P", P, "RH", 1)
            MW_ha = x_w_gas * MW_W + (1 - x_w_gas) * MW_A
            return {
                "T_DB (K)": T,
                "P (Pa)": P,
                "x_w (mol/mol)": Y,
                "x_w_gas (mol/mol)": x_w_gas,
                "h_gas (J/mol)": HAPropsSI("Hha", "T", T, "P", P, "RH", 1) * MW_ha,
                "s_gas (J/mol/K)": HAPropsSI("Sha", "T", T, "P", P, "RH", 1) * MW_ha,
                "v_gas (m3/mol)": HAPropsSI("Vha", "T", T, "P", P, "RH", 1) * MW_ha,
                "x_w_liq (mol/mol)": 1.0,
                "h_liq (J/mol)": PropsSI("H", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "s_liq (J/mol/K)": PropsSI("S", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "v_liq (m3/mol)": 1 / PropsSI("D", "T", T, "P", P, "water") * MW_W if T > Ttriple else 0,
                "x_gas (mol/mol)": (1 - Y) / (1 - x_w_gas),
                "RH": 1.0,
                "T_WB (K)": T
            }
    except Exception:
        return None

# --- Input Generator ---
def generate_random_inputs(size, T_bounds, P_bounds, Y_bounds):
    """
    Generates a list of random (T, P, Y) input tuples within specified bounds.

    Parameters:
        size (int): Number of samples
        T_bounds (tuple): Min and max temperature in Celsius
        P_bounds (tuple): Min and max pressure in kPa
        Y_bounds (tuple): Min and max mole fraction

    Returns:
        List of tuples (T [K], P [Pa], Y [mol/mol])
    """
    T_min, T_max = T_bounds
    P_min, P_max = P_bounds
    Y_min, Y_max = Y_bounds

    rand_T = np.random.uniform(T_min + 273.15, T_max + 273.15, size).round(6)
    rand_P = np.random.uniform(P_min * 1000, P_max * 1000, size).round(4)
    rand_Y = np.random.uniform(Y_min, Y_max, size).round(8)

    return list(zip(rand_T, rand_P, rand_Y))

# --- Internal Wrapper ---
def _get_data_point(args):
    """
    Internal helper that wraps compute_humid_air_properties and returns a DataFrame.
    Used for parallel execution.

    Parameters:
        args (tuple): A tuple of (T, P, Y)

    Returns:
        pd.DataFrame with one row or None if invalid
    """
    T, P, Y = args
    props = compute_humid_air_properties(T, P, Y)
    if props is None:
        return None
    return pd.DataFrame([props])

# --- Dataset Generator ---
def generate_HAprops_dataset(
    project_name: str = 'pysmo_humid_air_props',
    size: int = 100,
    T_bounds = (-30, 300), # degC
    P_bounds = (50, 200),  # kPa
    Y_bounds = (0.0, 0.2), # mol/mol
    save_csv: bool = True,
    serial: bool = False
) -> pd.DataFrame:
    """
    Generates a dataset of humid air properties using randomized inputs.

    Parameters:
        project_name (str): Output CSV filename prefix
        size (int): Number of data points to generate
        T_bounds (tuple): Temperature bounds in Celsius
        P_bounds (tuple): Pressure bounds in kPa
        Y_bounds (tuple): Water mole fraction bounds
        save_csv (bool): Whether to save the results to a CSV file
        serial (bool): Use serial (True) or parallel (False) processing

    Returns:
        pd.DataFrame containing valid humid air property data
    """
    inputs = generate_random_inputs(size, T_bounds, P_bounds, Y_bounds)

    startTime = timer()

    if serial:
        results = [_get_data_point(args) for args in inputs]
    else:
        with Pool() as p:
            results = p.map(_get_data_point, inputs)

    endTime = timer()

    df_list = [r for r in results if r is not None]
    df_final = pd.concat(df_list, ignore_index=True)

    print("Time taken: {:.2f}s".format(endTime - startTime))
    print("Valid points:", df_final.shape[0])

    if save_csv:
        df_final.to_csv(f"results/{project_name}_data.csv", index=False)

    return df_final

