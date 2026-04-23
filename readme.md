--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2026
--------------------------------------------------------------------------------


# PROJECT UPDATE:Too Many Cells processing Pipeline

## 1. OVERVIEW
This program aim's  to utilize the [too-many-cells](https://gregoryschwartz.github.io/too-many-cells/) spectral clustring
algorithm on [ImmuneDB](https://docs.immunedb.com/en/latest/) format BCR immune reportiores. Our [processing pipeline](https://github.com/dani-f99/tmc_preprocessing)
transforms the BCR genetic information into a set of unique 20-mer sequqnces clusterd by their 3-mer utilization using the spectral clustring algorithm from [too-many-cells](https://gregoryschwartz.github.io/too-many-cells/). In this way we identify common motiffs in the protein landscape of different B cell sub populations.



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
- [git](https://git-scm.com/)


## 3. USAGE GUIDE
1. Run the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) program (see it's documantation).
2. Clone the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing)() github reportiory (`git clone https://github.com/dani-f99/tmc_processing.git`).
3. Configure the `config.json` file (see section 5 of this readme).
4. See `tutorial_notebook.ipynb`, and run the first cell - which will automaticly create the folder srtucture.
5. Move the relevent files from the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) output to the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) input folders (for each analysis, will be expanded sepertly).
6. Run the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing), can use the `tutorial_notebook.ipynb`.


#### 3.1. CLUSTER ANALYSIS
For detailed examples please consult the `tutorial_notebook.ipynb`.
1. Move the `cluster_tree.json`, `clusters.csv` and `{x}-labels.csv` (x stands for metadata labels) files from the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) `tmc_output` folder to the appropriate `/input/{database-subject}` folder in the [tmc_processing](https://github.com/dani-f99/tmc_processing) cloned folder.
2. Initiate the clusters class by using the `source.clusters_analysis` class.
3. Get information about specific clusters of BCR cells by using the class method `source.clusters_analysis.cluster_analysis`.
4. Visuzlize the claster metadata fraction by using the command `source.clusters_analysis.lookup_node`.


#### 3.2. CLUMPINESS ANALYSIS
For detailed examples please consult the `tutorial_notebook.ipynb`.
1. Move the `clusters.csv` and `labels_joined.csv` files from the [tmc_preprocessing](https://github.com/dani-f99/tmc_preprocessing) `tmc_output` folder to the appropriate `/input/{database-subject}` folder in the [tmc_processing](https://github.com/dani-f99/tmc_processing) cloned folder.
2. Use the function `source.clumpiness_analysis.clumpiness_json` to create appropriate [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) json input.
3. Install [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) (for detailed installation tutorial see section 4).
4. Use the command: `cat input/subjX/find_clumpiness_input.json | find-clumpiness --format JSON > output/subjX/clumpiness_data.csv` to create clumpiness analysis based on the input data, this may require high amounts of RAM (path in the command is an example, can modify as you see fit).
5. Move the [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) `clumpiness_data.csv` output into the [tmc_processing](https://github.com/dani-f99/tmc_processing) `input` folder.
6. Run the function `source.clumpiness_analysis.clumpiness_heatmap` to plot clumpiness heatmap.


#### 3.3. CLONES ANALYSIS




## 4. FIND CLUMPINESS INSTALLATION
[find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) runs nativly on linux machine, if windows is insalled use the compatability layer `wsl2`.
1. clone the [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) into your machine by using the command `git clone https://github.com/GregorySchwartz/find-clumpiness.git`.
2. In the [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) folder, 
3. Install haskell by using the command `apt install haskell-stack` (may require sudo).
4. In the [find-clumpiness](https://github.com/GregorySchwartz/find-clumpiness) folder, run the command `stack install find-clumpiness`.
5. if `Warning: Installation path /home/daniel/.local/bin not found on the PATH environment variable.` appears run the following commands: `echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc` and `source ~/.bashrc`.
6. Verify installation by running the command `which find-clumpiness` (expected output: `../.local/bin/find-clumpiness`).


## 5. CONFIG.JSON STRUCTURE AND CONFIGURATION
The purpose of the cofig file is to make it easier and faster to utilize this program, we define two variables 
which are database and subjects:
- `database` - database name (from the previous analysis - see 3.Usage guide).
- `subjects` - subjects to be analyzed. need to be a string of numbers with `,` delimiter (for example: "3,4,5").



## 6. DIRECTORY STRUCTURE
The program uses the following folder structure, where in each main folder there is a database-subject specific sub-folder
with it's dedicated files.
 - `input` -> Input folder. Used to store the nessesery files used by this program.
 - `output` -> Output folder. Used to store the output figures and datasets of this program.
 
	
	
## 7. RESOURCES
- tmc_preprocessing: https://github.com/dani-f99/tmc_preprocessing
- too-many-cells-python: https://gregoryschwartz.github.io/too-many-cells/
- find-clumpiness: https://github.com/GregorySchwartz/find-clumpiness
- ImmuneDB: https://docs.immunedb.com/en/latest/