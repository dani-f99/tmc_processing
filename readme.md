--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2026
--------------------------------------------------------------------------------


# PROJECT UPDATE:Too Many Cells processing Pipeline

## 1. OVERVIEW
This program aim is to utilize the [too-many-cells](https://gregoryschwartz.github.io/too-many-cells/) spectral clustring
algorithm on [ImmuneDB](https://docs.immunedb.com/en/latest/) format BCR immune reportiores. Our [processing pipeline](https://github.com/dani-f99/tmc_preprocessing)
prepares the BCR genetic information as unique 20-mer sequqnces crossed with 3-mer utilization whitin them in order to unveil genetic motifs common
in different BCR sub-populations by the spectral clustring algorithm. The [too-many-cells](https://gregoryschwartz.github.io/too-many-cells/) output is used as input for this program.


## 2. PREREQUISITES
This program works on ImmuneDB databse (see resources section)

Python Libraries:
- `Pandas`
- `NumpPy`
- `Matplotlib`
- `Seaborn`

Programs:
- [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing)
- [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness)


## 3. USAGE GUIDE
1. Run the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) program (see it's documantation).
2. Clone the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing)() github reportiory (`git clone https://github.com/dani-f99/tmc_processing.git`).
3. Configure the `config.json` file (see section 4 of this readme).
4. See `tutorial_notebook.ipynb`, and run the first cell - which will automaticly create the folder srtucture.
5. Move the relevent files from the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) output to the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing)() input folders (for each analysis, will be expanded sepertly).
6. Run the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing), can use the `tutorial_notebook.ipynb`.


#### 3.1. CLUSTER ANALYSIS


#### 3.2. CLUMPINESS ANALYSIS


#### 3.3. CLONES ANALYSIS



## 4. CONFIG.JSON STRUCTURE AND CONFIGURATION
The purpose of the cofig file is to make it easier and faster to utilize this program, we define two variables 
which are database and subjects:
- `database` - database name (from the previous analysis - see 3.Usage guide).
- `subjects` - subjects to be analyzed. need to be a string of numbers with `,` delimiter (for example: "3,4,5").



## 5. DIRECTORY STRUCTURE
The program uses the following folder structure, where in each main folder there is a database-subject specific sub-folder
with it's dedicated files.
 - `input` -> Input folder. Used to store the nessesery files used by this program.
 - `output` -> Output folder. Used to store the output figures and datasets of this program.
 
	
	
## 6. RESOURCES
- tmc_preprocessing: https://github.com/dani-f99/tmc_preprocessing
- too-many-cells-python: https://gregoryschwartz.github.io/too-many-cells/
- find-clumpiness: https://github.com/GregorySchwartz/find-clumpiness
- ImmuneDB: https://docs.immunedb.com/en/latest/