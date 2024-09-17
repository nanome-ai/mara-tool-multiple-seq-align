import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics
import biotite.sequence.io.fasta as fasta
import numpy as np
from datetime import datetime
import pickle

def create_timestamp() -> str:
  # helper function to create a unique timestamp
  dt = str(datetime.now())
  return dt.replace("-", "_").replace(":", "_").replace(" ", "_")

def run(seq_fasta: str, reorder=False):
    with open(seq_fasta, "r") as file:
        fasta_file = fasta.FastaFile.read(file)

    # Create Sequence objects for the alignment
    seqs = []
    for header, sequence in fasta_file.items():
        seqs.append(seq.ProteinSequence(sequence))  

    # make sure there are more than 1 sequences
    assert len(seqs) > 1, "Error! The input fasta file does not contain multiple sequences!"

    # Actual multiple sequence alignment step
    alignment, order, guide_tree, distance_matrix = align.align_multiple(
        seqs,
        matrix=align.SubstitutionMatrix.std_protein_matrix(),
        gap_penalty=-5,
    )
    # Order alignment according to guide tree if needed and draw
    labels = np.array(list(fasta_file.keys()))
    if reorder:
        labels = labels[order]
        alignment = alignment[:, order.tolist()]

    # fig, ax = plt.subplots(figsize=(6.0, 3.0), constrained_layout=True)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_axes([0.2, 0.1, 0.7, 0.8])
    ax = plt.gca()
    graphics.plot_alignment_similarity_based(
        ax, alignment, 
        symbols_per_line=min(50, len(alignment)), show_numbers=True, symbol_size=10,
        labels=labels, label_size=12,
        show_line_position=True
    )

    # Save the figure as a svg file, and the alignment object as pkl for later use
    output_name = f"alignment.{create_timestamp()}"

    plt.savefig(output_name + ".svg")
    print(f'Alignment image saved as {output_name}.svg')

    with open(output_name + ".fasta", "w") as file:
        for i in range(len(alignment.sequences)):
            file.write(f">{labels[i]}\n{alignment._gapped_str(i)}\n")
    print("Aligned sequences in fasta format are saved as:", output_name + ".fasta")

    with open(output_name + ".pkl", "wb") as file:
        pickle.dump((alignment, labels, guide_tree.to_newick()), file)

    print(f'Intermediate results saved to {output_name}.pkl. Please update this file in context for further analysis.')

    # Close the figure to free up memory
    plt.close()

   

if __name__ == '__main__':
    run("example_data/two_seq.fasta")