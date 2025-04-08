from compute_HAprops import generate_HAprops_dataset
from train_surrogate import train_HA_pp


if __name__ == "__main__":
    project_name='pysmo_humid_air_props'
    generate_HAprops_dataset(project_name, size=500000)
    train_HA_pp(project_name, n_inputs=3, sample_rows=4000, model_type='rbf')