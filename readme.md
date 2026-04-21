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
#### 3.1. CLUSTER ANALYSIS


#### 3.2. CLUMPINESS ANALYSIS


#### 3.3. CLONES ANALYSIS


## 4. CONFIG.JSON STRUCTURE AND CONFIGURATION



## 5. DIRECTORY STRUCTURE
The program uses the following folder structure: 

	
	
## 6. RESOURCES
- tmc_preprocessing: https://github.com/dani-f99/tmc_preprocessing
- too-many-cells-python: https://gregoryschwartz.github.io/too-many-cells/
- find-clumpiness: https://github.com/GregorySchwartz/find-clumpiness
- ImmuneDB: https://docs.immunedb.com/en/latest/
