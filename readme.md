Toos:

* Search similar sequence using NCBI blast
* Combine multiple fasta files into one
* Multiple sequence alignment and visualization
* Visualize part of the aligned sequences
* Cluster aligned sequences based on similarity
* Analyze differences of two aligned sequences



# Multiple sequence alignment and visualization (Informatics)

Take multiple sequences and align them with visual output. 

## Input
seq_fasta - file path to a single fasta file containing multiple sequences to be aligned
reorder - whether the sequences should be reordered based on the similarity after alignment. Default to False

## Output
* An image showing the aligned sequences, which are highlighted by amino acid similarity.
* Aligned sequences in the fasta format
* A pickle file containing the alignment object, the labels and the newick representation of the alignment tree, which can be used by other tools

**Important! Please display the image in the response!**


# Cluster aligned sequences based on similarity (Informatics)

Using an intermediate result from a multiple sequence alignment, cluster the sequences based on their similarity and generate a dendrogram.

## Input
intermediate_pkl_file - file path to the pickle file created from the "Multiple sequence alignment and visualization" tool

## Output
* A dendrogram showing the clustering of the sequences 

# Analyze differences of two aligned sequences (Informatics)

Give a written summary of the differences between two aligned sequences. In the first part, sequences are aligned and gaps are indicated by dashes. Different residues are represented by lower-case letters. In the second part, all different residues are enumerated. You should use the full names of the amino acids to explain the differences.

## Input
intermediate_pkl_file - file path to the pickle file created from the "Multiple sequence alignment and visualization" tool. The alignment object should contain exactly two sequences.

## Output
The differences in the amino acid sequences between the two aligned sequences.

# Setup
### requirements
biotite
matplotlib
numpy