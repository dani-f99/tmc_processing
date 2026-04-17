# Import required packages
import importlib.metadata
import pandas as pd
import json
import os


##############################################################
# Custom function that cheeck if python packages are installed
def check_packages(package_names):
    status = {}

    n_total = len(package_names)
    n_installed = 0
    package_df = pd.DataFrame(index=package_names, columns=["installed","not_installed"], data=0)

    for pkg in package_names:
        try:
            # metadata.version returns the version string if installed
            dist_version = importlib.metadata.version(pkg)
            status[pkg] = f"Installed (v{dist_version})"
            package_df.loc[pkg, "installed"] = 1
            n_installed += 1

        except importlib.metadata.PackageNotFoundError:
            status[pkg] = "Not Installed !!!"
            package_df.loc[pkg, "not_installed"] = 1
    
    print(f"> {n_installed}/{n_total} packages are installed:")

    for i in status:
        print(f"{i} package is {status[i]}")

    return package_df


################################################################################$$##########
# Reading information from json file. Used to extract the parameters from the `config.json`.
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open(path) as config:
        config_f = json.load(config)

    return config_f


#####################################################
# Creating folder according to the and program scheme
def create_folders(req_folders : list = ["temp_data", "tms_input", "reports"]):
    """
    req_folders : str -> required folders path, if subfolder exsits input '\\' between folders.
    """

    for folder in req_folders:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
            print(f"> folder `{folder}` was created.")

        else:
            print(f"> folder `{folder}` exists, continuing.")
