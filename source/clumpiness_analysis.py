####################################
from source.helpers import read_json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
import sys
import os


################################################
def clumpiness_json(path_node_info : str = None,
                    path_labels : str = None):
    """
    path_node_info : str -> alternative path of the `cells_clusters_info.csv` file, defualt path is in `intpu/db-subject/cells_clusters_info.csv`.
    path_labels : str -> alternative path for the `labels.csv` file, defualt path is in `intpu/db-subject/labels.json`.
    """

    # Config information
    config = read_json()
    database , subject = config["database"], config["subjects"]

    # Increase recursion limit for the deep tree analysis
    sys.setrecursionlimit(20000)

    # 1. Load the data
    dfs_dict = {}
    for path, label in zip([path_node_info, path_labels],["cells_clusters_info.csv", "labels.csv"]):
        if path is None:
            dfs_dict[label] = pd.read_csv(os.path.join("input", f"{database}-subject{subject}", f"{label}"))
        else:
            dfs_dict[label] = pd.read_csv(path)

    clusters_df = dfs_dict["cells_clusters_info.csv"]
    labels_df = dfs_dict["labels.csv"]

    # Use the first column (barcodes) and rename for clarity
    clusters_df.rename(columns={clusters_df.columns[0]: 'item'}, inplace=True)

    # 2. Map items to their specific labels
    item_label_map = dict(zip(labels_df['item'], labels_df['label']))

    # 3. Build adjacency maps for the sp tree
    cluster_to_children = {} # parent cluster -> set of child clusters
    cluster_to_items = {}    # leaf cluster -> list of (barcode, label)

    for _, row in clusters_df.iterrows():
        path = str(row['sp_path']).split('/') # e.g., ["25523", ..., "0"]
        item = row['item']
        label = item_label_map.get(item)
        
        # Map the item to its most specific sp cluster (the first in the path)
        leaf_cluster = path[0]
        if leaf_cluster not in cluster_to_items:
            cluster_to_items[leaf_cluster] = []
        cluster_to_items[leaf_cluster].append((item, label))
        
        # Build hierarchy: path[i] is a child of path[i+1]
        for i in range(len(path) - 1):
            child = path[i]
            parent = path[i+1]
            if parent not in cluster_to_children:
                cluster_to_children[parent] = set()
            cluster_to_children[parent].add(child)

    # 4. Recursive function to build the [metadata, children] tuple format
    def build_clumpiness_node(node_id):
        # The first element is the metadata object
        metadata = {
            "nodeID": str(node_id),
            "nodeLabels": [] # Internal nodes usually have empty labels
        }
        
        # The second element is the list of children nodes
        children = []
        
        # Add child clusters
        if node_id in cluster_to_children:
            for child_id in cluster_to_children[node_id]:
                children.append(build_clumpiness_node(child_id))
                
        # Add individual cell items as leaf nodes
        if node_id in cluster_to_items:
            for item_id, label in cluster_to_items[node_id]:
                leaf_metadata = {
                    "nodeID": str(item_id),
                    "nodeLabels": [str(label)] if label else []
                }
                # A leaf is a 2-element list with an empty children list: [meta, []]
                children.append([leaf_metadata, []])
                
        return [metadata, children]

    # 5. Generate tree from root '0' and save
    final_tree = build_clumpiness_node('0')

    output_path = os.path.join("output", f"{database}-subject{subject}")
    with open(os.path.join(output_path, 'find_clumpiness_input.json'), 'w') as f:
        json.dump(final_tree, f)

    print(f"`find_clumpiness_input.json` Saved at {output_path}.")


#############################################
def clumpiness_heatmap(dataset : str = None,
                       plot_name : str = None):
    
    """
    dataset : str -> csv file location, if None then will loog for `clumpiness_data.csv` at the input folder.
    plot_name : str -> descriptive name for the clumpiness heatmap plot file.
    """

    db, subj = read_json()["database"], read_json()["subjects"]
    plot_path = os.path.join("output", f"{db}-subject{subj}", f"clumpiness_heatmap-{plot_name}.png")

    if dataset is not None:
        # Import clumpiness from csv.
        try:
            dataset = pd.read_csv(dataset)
        except:
            raise Exception("Invalid clumpiness csv path.")

    # Import clumpiness csv from file at the 'data' folder.
    else:
        data_file = os.path.join("input", f"{db}-subject{subj}", "clumpiness_data.csv")
        dataset = pd.read_csv(data_file)

    # Getting axis names (mirrored x and y for heatmap).
    axis_names = np.sort(np.unique(np.concatenate((dataset.iloc[:,0].values, 
                                                dataset.iloc[:,1].values))))
    
    # Creating template dataframe.
    clumpiness_df = pd.DataFrame(index = axis_names, columns=axis_names, data=np.nan)
    
    # Assigning values from the sparce matrix to the template matrix.
    label_col0 , label_col1 =  dataset.columns[:-1]
    
    for i in axis_names:
        for j in axis_names:
            try:
                cell_data = dataset.loc[(dataset[label_col0] == i) & (dataset[label_col1] == j), "value"].values[0]
            except:
                cell_data = np.nan

            clumpiness_df.loc[i, j] = cell_data

    sns.heatmap(clumpiness_df)
    plt.savefig(plot_path, bbox_inches='tight')
    print(f"Plot {plot_name} saved to `{plot_path}`.")
    plt.show()