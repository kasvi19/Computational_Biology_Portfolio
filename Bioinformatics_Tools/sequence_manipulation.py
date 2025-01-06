from collections import defaultdict

# Complement dictionary that defaults to the input nucleotide if not found
complement = defaultdict(lambda x: x, {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
})

def reverse_complement(seq, complement):
    reverse_complement_seq = "".join(complement.get(nt, nt) for nt in reversed(seq))
    return reverse_complement_seq






