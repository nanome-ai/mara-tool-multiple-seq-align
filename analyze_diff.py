import pickle

aa_map = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic acid',
    'C': 'Cysteine',
    'Q': 'Glutamine',
    'E': 'Glutamic acid',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'V': 'Valine',
    'B': 'Aspartic Acid or Asparagine (ambiguous)',
    'Z': 'Glutamic Acid or Glutamine (ambiguous)',
    'J': 'Leucine or Isoleucine (ambiguous)',
    'X': 'Unknown or any amino acid'
}

def run(intermediate_pkl_file):
    with open(intermediate_pkl_file, "rb") as file:
        alignment, labels, alignment_tree_newick = pickle.load(file)

    if len(alignment.sequences) != 2:
        raise ValueError("This tool can be only used for pairwise sequence alignment! \
Please make sure the input to the alignment are two sequences.")
    

    # Print the differences between the sequences
    print(f'Differences between the sequences:\n')
    trace = list(map(tuple, alignment.trace))
    aligned_seq1 = ''.join(alignment.sequences[0][i] if i != -1 else '-' for i in [t[0] for t in trace])
    aligned_seq2 = ''.join(alignment.sequences[1][i] if i != -1 else '-' for i in [t[1] for t in trace])
    diff_seq = ""
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b:
            diff_seq += "-"
        else:
            diff_seq += "|"
    
    # Split the sequences into chunks of 80 characters for better readability
    chunk_size = 50
    aligned_seq1_chunks = [aligned_seq1[i:i+chunk_size] for i in range(0, len(aligned_seq1), chunk_size)]
    aligned_seq2_chunks = [aligned_seq2[i:i+chunk_size] for i in range(0, len(aligned_seq2), chunk_size)]
    diff_seq_chunks = [diff_seq[i:i+chunk_size] for i in range(0, len(diff_seq), chunk_size)]
    
    # Collect differences
    differences = []
    
    # Print the aligned sequences and differences in chunks
    for i, (seq1_chunk, seq2_chunk) in enumerate(zip(aligned_seq1_chunks, aligned_seq2_chunks), start=1):
        start_pos = (i - 1) * chunk_size + 1
        end_pos = min(i * chunk_size, len(aligned_seq1))
        print(f"Sequence segment {i} (positions {start_pos}-{end_pos})\n")
        
        # Generate the modified sequence chunks with uppercase and lowercase letters
        seq1_modified = ''
        seq2_modified = ''
        for j, (a, b) in enumerate(zip(seq1_chunk, seq2_chunk), start=start_pos):
            if a == b:
                if a == '-':
                    seq1_modified += '-'
                    seq2_modified += '-'
                else:
                    seq1_modified += a.upper()
                    seq2_modified += b.upper()
            else:
                if a == '-' or b == '-':
                    seq1_modified += a
                    seq2_modified += b
                else:
                    seq1_modified += a.lower()
                    seq2_modified += b.lower()
                    differences.append((j, a, b))
        
        print(f'{seq1_modified}\n')
        print(f'{seq2_modified}\n')
    
    # Print the collected differences
    if differences:
        print("Differences in chains found at position (aligned numbering):")
        print(f"|Position|{labels[0]}|{labels[1]}|\n|---|---|---|")
        for pos, aa1, aa2 in differences:
            print(f"|{pos}|{aa_map[aa1]} ({aa1})|{aa_map[aa2]} ({aa2})|")
    else:
        print("No differences found between the chains.")


if __name__ == '__main__':
    run("example_data/two_seq.pkl")