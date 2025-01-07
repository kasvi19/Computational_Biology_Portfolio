from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
import csv
import os

# Complement dictionary that defaults to the input nucleotide if not found
complement = defaultdict(lambda x: x, {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
})

# Function to find reverse complement of a sequence
def reverse_complement(records, output_fasta):
    reverse_complemented_records = []
    for record in records:
        seq = str(record.seq)
        reverse_complement_seq = "".join(complement.get(nt, nt) for nt in reversed(seq.upper()))
        reverse_complemented_records.append(SeqRecord(Seq(reverse_complement_seq), id=record.id, description="Reverse Complement"))
    
    # Write reverse complement sequences to FASTA file
    with open(output_fasta, "w") as rc_fasta_file:
        SeqIO.write(reverse_complemented_records, rc_fasta_file, "fasta")
    return reverse_complement_seq


# Function to find mRNA sequence of strand, treats it as coding strand by default
def mrna(records, output_fasta, is_template_strand=False):
    mrna_records = []
    for record in records:
        seq = str(record.seq)
        if is_template_strand:
            seq = reverse_complement(seq, complement)
        mrna_seq = seq.replace('T','U')
        mrna_records.append(SeqRecord(Seq(mrna_seq), id=record.id, description="mRNA"))
    
    # Write mRNA sequences to FASTA file
    with open(output_fasta, "w") as mrna_fasta_file:
        SeqIO.write(mrna_records, mrna_fasta_file, "fasta")
    return mrna_seq


# Function to translate DNA sequences to protein sequences and calculate stats like stop codons, total amino acids, invalid codons
def translation(input_data, output_fasta, output_csv):
    translated_records = []                                  # List to store translated sequences
    standard_table = CodonTable.unambiguous_dna_by_id[1]

    # Parse input_data into records
    if isinstance(input_data, str):                          # If a file path is given
        records = SeqIO.parse(input_data, "fasta")
    else:                                                    # Assume it's a list of SeqRecord objects
        records = input_data

    # Initialize stats dictionary
    stats = {
        "total_sequences": 0,
        "total_stop_codons": 0,
        "total_aa": 0,
        "invalid_codons": 0,
    }

    # Process each sequence individually
    for record in records:
        seq = str(record.seq).upper()
        # Any trailing nucleotides that do not form a codon are removed
        seq = seq[:-(len(seq) % 3)] if len(seq) % 3 != 0 else seq

        stats["total_sequences"] += 1
        translated_sequence = ""
        try:
            # Translate the sequence including stop codons
            protein = record.seq.translate(to_stop=False, table=standard_table)
            translated_sequence = str(protein)
            stats["total_stop_codons"] += translated_sequence.count("*")
            stats["total_aa"] += len(protein)
        except Exception as e:
            print(f"Error translating sequence {record.id}: {e}")
            stats["invalid_codons"] += len(seq) // 3                              # Assuming all codons failed
        translated_records.append(SeqRecord(Seq(translated_sequence), id=record.id, description=""))

    # Translated sequences to output FASTA file
    with open(output_fasta, "w") as fasta_file:
        SeqIO.write(translated_records, fasta_file, "fasta")
    # CSV file for stats
    with open(output_csv, mode="w", newline="") as csvfile:
        fieldnames = ['Total Sequences', 'Total Stop Codons', 'Total AA', 'Invalid Codons']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({
            'Total Sequences': stats["total_sequences"],
            'Total Stop Codons': stats["total_stop_codons"],
            'Total AA': stats["total_aa"],
            'Invalid Codons': stats["invalid_codons"]
        })

    print(f"Translation complete. Stats and sequences written to {output_fasta} and {output_csv}.")
    return translated_records, stats

# Wrapper function to call all the sequence manipulation functions
def sequence_manipulation(input_file, output_fasta, output_csv):
    records = list(SeqIO.parse(input_file, "fasta"))
    # Output filenames with appropriate prefixes
    reverse_complement_output = os.path.join(os.path.dirname(output_fasta), f"rc_{os.path.basename(output_fasta)}")
    mrna_output = os.path.join(os.path.dirname(output_fasta), f"mrna_{os.path.basename(output_fasta)}")
    translated_output = os.path.join(os.path.dirname(output_fasta), f"translated_{os.path.basename(output_fasta)}")
    reverse_complement(records, reverse_complement_output)
    mrna(records, mrna_output, is_template_strand=False)
    translation(input_file, translated_output, output_csv)
    
    print(f"Reverse complement, mRNA sequence, and translation complete. Results saved to {reverse_complement_output}, {mrna_output}, and {translated_output}.")

# Usage
input_file = os.path.join(os.path.dirname(__file__), '../Sample_Data/coronavirus.fasta')
output_fasta = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/coronavirus_seq_manipulation.fasta')
output_csv = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/coronavirus_translation_stats.csv')
sequence_manipulation(input_file, output_fasta, output_csv)

