import pandas as pd
import matplotlib
matplotlib.use('Agg')

from typing import Tuple, List, Union, Optional
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import (
    PysmoPolyTrainer, PysmoRBFTrainer, PysmoKrigingTrainer, PysmoSurrogate
)
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D, surrogate_parity, surrogate_residual,
)

# ----------------------------------------------------------------------------------
# Constants and header mapping

HEADER_MAP = {
    'T_DB (K)': 'temperature',
    'P (Pa)': 'pressure',
    'x_w (mol/mol)': 'mole_frac_1',
    'x_w_gas (mol/mol)': 'mole_frac_1_gas_phase',
    'h_gas (J/mol)': 'enth_mol_gas_phase',
    's_gas (J/mol/K)': 'entr_mol_gas_phase',
    'v_gas (m3/mol)': 'vol_mol_gas_phase',
    'x_w_liq (mol/mol)': 'mole_frac_1_liq_phase',
    'h_liq (J/mol)': 'enth_mol_liq_phase',
    's_liq (J/mol/K)': 'entr_mol_liq_phase',
    'v_liq (m3/mol)': 'vol_mol_liq_phase',
    'x_gas (mol/mol)': 'mol_frac_gas_phase',
    'RH': 'relative_humidity',
    'T_WB (K)': 'temperature_wet_bulb',
}

# surr_property_map = {
#     'temperature': 'K',
#     'pressure': 'Pa',
#     'mole_frac_1': 'mol/mol',
#     'mole_frac_1_gas_phase': 'mol/mol',
#     'enth_mol_gas_phase': 'J/mol',
#     'entr_mol_gas_phase': 'J/mol/K',
#     'vol_mol_gas_phase': 'm3/mol',
#     'mole_frac_1_liq_phase': 'mol/mol',
#     'enth_mol_liq_phase': 'J/mol',
#     'entr_mol_liq_phase': 'J/mol/K',
#     'vol_mol_liq_phase': 'm3/mol',
#     'mol_frac_gas_phase': 'mol/mol',
#     'relative_humidity': 'mol/mol',
#     'temperature_wet_bulb': 'K',
# }

COLS = len(HEADER_MAP)

# ----------------------------------------------------------------------------------
# Data utilities

def get_training_data_from_file(filename: str, sample_rows: int = 1000) -> pd.DataFrame:
    df = pd.read_csv(filename, header=0, index_col=False)
    if df.empty:
        raise ValueError("CSV file is empty.")
    for col in df.columns:
        if col in HEADER_MAP:
            df.rename(columns={col: HEADER_MAP[col]}, inplace=True)
    return df.sample(n=min(sample_rows, len(df)))

# ----------------------------------------------------------------------------------
# Plotting

def get_model_validation_plots(
    poly_surr: PysmoSurrogate,
    data: pd.DataFrame,
    model_name: str,
    data_tag: str = 'val'
) -> None:
    if data.empty:
        raise ValueError("Cannot generate plots for empty DataFrame.")
    surrogate_scatter2D(poly_surr, data, filename=f"results/{model_name}_{data_tag}_scatter2D.pdf")
    surrogate_parity(poly_surr, data, filename=f"results/{model_name}_{data_tag}_parity.pdf")
    surrogate_residual(poly_surr, data, filename=f"results/{model_name}_{data_tag}_residual.pdf")

# ----------------------------------------------------------------------------------
# Surrogate training

def _get_trainer(
    model_type: str,
    input_labels: List[str],
    output_labels: List[str],
    training_data: pd.DataFrame
) -> Union[PysmoPolyTrainer, PysmoRBFTrainer, PysmoKrigingTrainer]:
    trainers = {
        'poly': PysmoPolyTrainer,
        'rbf': PysmoRBFTrainer,
        'kriging': PysmoKrigingTrainer
    }
    trainer_cls = trainers.get(model_type.lower())
    if trainer_cls is None:
        raise ValueError(f"Unsupported model_type '{model_type}'. Choose from {set(trainers)}.")

    trainer = trainer_cls(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe=training_data
    )

    if model_type == 'poly':
        trainer.config.extra_features = [
            "pressure*temperature*temperature",
            "pressure*pressure*temperature*temperature",
            "pressure*pressure*temperature",
            "pressure/temperature",
            "temperature/pressure",
        ]
        trainer.config.maximum_polynomial_order = 3 if training_data.shape[0] < 100 else 5
        trainer.config.multinomials = True

    if hasattr(trainer.config, "training_split"):
        trainer.config.training_split = 0.8
    if hasattr(trainer.config, "number_of_crossvalidations"):
        trainer.config.number_of_crossvalidations = 10

    return trainer


# ----------------------------------------------------------------------------------
# Main API

def train_surrogate_model(
    df: pd.DataFrame,
    n_inputs: int,
    xmin: Tuple[float, ...],
    xmax: Tuple[float, ...],
    model_type: str = 'rbf',
    output_filename: str = 'surrogate_model',
    config: Optional[dict] = None
) -> Tuple[PysmoSurrogate, pd.DataFrame, pd.DataFrame]:

    if df.empty:
        raise ValueError("Training DataFrame is empty.")

    input_labels = list(df.columns)[0:n_inputs]
    output_labels = list(df.columns)[n_inputs:COLS]

    if len(xmin) != len(input_labels) or len(xmax) != len(input_labels):
        raise ValueError("Length of xmin and xmax must match number of input features.")

    train_df, val_df = split_training_validation(df, 0.5, seed=len(df))
    trainer = _get_trainer(model_type, input_labels, output_labels, train_df)

    if config:
        for key, val in config.items():
            setattr(trainer.config, key, val)

    trained_model = trainer.train_surrogate()
    input_bounds = {name: (xmin[i], xmax[i]) for i, name in enumerate(input_labels)}
    surrogate = PysmoSurrogate(trained_model, input_labels, output_labels, input_bounds)

    surrogate.save_to_file(f"results/{output_filename}_surrogate.json", overwrite=True)

    return surrogate, train_df, val_df

def train_HA_pp(project_name: str = "pysmo_humid_air", n_inputs: int = 3, sample_rows = 500, model_type = 'rbf') -> None:
    """
    Example script execution for training a surrogate model on humid air data.

    Loads data, scales it, trains the selected surrogate model, and saves it.
    Also generates training and validation plots.

    Change `model_type` to 'poly', 'rbf', or 'kriging' as needed.
    """
    
    df = get_training_data_from_file(f"results/{project_name}_data.csv", sample_rows)
    # df = df[['temperature', 'pressure', 'mole_frac_1', 'enth_mol_gas_phase']]

    # Define physical bounds
    xmin = [273.15 - 70, 50000, 0]
    xmax = [273.15 + 400, 200000, 0.5]

    # Train and export model
    surrogate, data_train, data_val = train_surrogate_model(
        df=df,
        n_inputs=n_inputs,
        xmin=xmin,
        xmax=xmax,
        model_type=model_type,
        output_filename=project_name
    )

    # Generate validation plots
    get_model_validation_plots(surrogate, data_train, project_name, 'train')
    get_model_validation_plots(surrogate, data_val, project_name, 'val')

    print("✅ Surrogate training and export complete.")

if __name__ == "__main__":
    train_HA_pp()