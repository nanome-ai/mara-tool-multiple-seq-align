import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics
import biotite.sequence.phylo as phylo
import numpy as np
from datetime import datetime
import pickle

def create_timestamp() -> str:
  # helper function to create a unique timestamp
  dt = str(datetime.now())
  return dt.replace("-", "_").replace(":", "_").replace(" ", "_")

def run(intermediate_pkl_file):
    with open(intermediate_pkl_file, "rb") as file:
        alignment, labels, alignment_tree_newick = pickle.load(file)

    alignment_tree = phylo.Tree.from_newick(alignment_tree_newick)

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_axes([0.06, 0.1, 0.7, 0.85])
    graphics.plot_dendrogram(
        ax, alignment_tree, 
        labels=labels,
        label_size=18,
        linewidth=2,
    )

    # Save the figure as a svg file, and the alignment object as pkl for later use
    output_name = f"dendrogram.{create_timestamp()}"

    plt.savefig(output_name + ".svg")
    print(f'Clustered dendrogram image saved to {output_name}.svg')

    plt.close()

   

if __name__ == '__main__':
    run("example_data/short.pkl")