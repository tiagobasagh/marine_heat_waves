import os

DIR_BASE = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]

DIR_SST = os.path.join(DIR_BASE, "SST")
DIR_NOTEBOOKS = os.path.join(DIR_BASE, "notebooks")
DIR_DATASETS = os.path.join(DIR_NOTEBOOKS,"datasets")

def save_dataset_as_nc(ds, dir_to_save=None, name_file="example.nc"):
    """
    
    """
    if not dir_to_save:
        dir_to_save = DIR_DATASETS

    ds.to_netcdf(os.path.join(dir_to_save, name_file))


def save_dataframe_as_csv(df, dir_to_save=None, name_file="exaple.csv"):
    """
    """
    if not dir_to_save:
        dir_to_save = DIR_DATASETS
    df.to_csv(join_path(dir_to_save, name_file), index=False)


def join_path(path1, path2):
    """
    
    """
    return os.path.join(path1, path2)


def exist_file(path):
    """
    
    """
    return os.path.isfile(path)