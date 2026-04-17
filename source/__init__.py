from source.helpers import check_packages, create_folders, read_json
import os

# Cheeking if needed
check_packages(["pandas", "numpy", "matplotlib", "seaborn"])

# Importing databse and subject information from the config.json
config = read_json()
database = config["database"]
subejcts = config["subjects"].split(",")

# Creating folders - if dont exists
create_folders(["input", "output"]) 
create_folders([os.path.join("input",f"{database}-subject{subject}") for subject in subejcts]) 
create_folders([os.path.join("output",f"{database}-subject{subject}") for subject in subejcts])

