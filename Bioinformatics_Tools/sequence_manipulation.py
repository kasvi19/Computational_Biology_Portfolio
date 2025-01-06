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
stats = {
    "total_sequences": 0,
    "total_codons": 0,
    "total_aa": 0,
    "total_stop_codons": 0,
    "invalid_codons": 0
}

def translation(input_file, output_fasta, output_csv):
    translated_records = []    # List to store translated sequences
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    
    # Write the stats to a CSV file
    with open(output_csv, mode="w", newline="") as csvfile:
        fieldnames = ['Total Sequences', 'Total Codons', 'Total AA', 'Stop Codons', 'Invalid Codons']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
    
        # Parse the fasta file and translate the sequences
        for record in SeqIO.parse(input_file, "fasta"):
            seq = str(record.seq).upper()
            stats["total_sequences"] += 1
            translated_sequence = ""
            stop_codons = 0
            invalid_codons = 0
            try:
                # Translate the sequence until the end, regardless of stop codons
                protein = record.seq.translate(to_stop=False)
                translated_sequence = str(protein)  # Convert the translated protein to string

                # Count the number of stop codons and invalid codons
                stop_codons = translated_sequence.count("*")
                stats["total_stop_codons"] += stop_codons
                stats["total_aa"] += len(translated_sequence)
                invalid_codons = translated_sequence.count("X")
                
            except Exception as e:
                print(f"Error translating sequence {record.id}: {e}")
                invalid_codons = len(seq) // 3
                continue
            translated_records.append(SeqRecord(Seq(translated_sequence), id=record.id, description=""))
        
        # Write translated sequences to output FASTA file
        with open(output_fasta, "w") as fasta_file:
            SeqIO.write(translated_records, fasta_file, "fasta")
        
        # Write stats to the CSV file
        stats_to_write = {
            'Total Sequences': stats["total_sequences"],
            'Total Codons': stats["total_codons"],  
            'Total AA': stats["total_aa"],
            'Stop Codons': stats["total_stop_codons"],
            'Invalid Codons': stats["invalid_codons"]
        }
        writer.writerow(stats_to_write)
        print(f"Translation complete. Stats and sequences written to {output_fasta} and {output_csv}.")

# Wrapper function to call all the sequence manipulation functions
def sequence_manipulation(input_file, output_fasta, output_csv):
    records = list(SeqIO.parse(input_file, "fasta"))
    
    # Define output filenames with appropriate prefixes
    reverse_complement_output = os.path.join(os.path.dirname(output_fasta), f"rc_{os.path.basename(output_fasta)}")
    mrna_output = os.path.join(os.path.dirname(output_fasta), f"mrna_{os.path.basename(output_fasta)}")
    translated_output = os.path.join(os.path.dirname(output_fasta), f"translated_{os.path.basename(output_fasta)}")
    
    # Call the individual functions with updated output filenames
    reverse_complement(records, reverse_complement_output)
    mrna(records, mrna_output, is_template_strand=False)
    translation(input_file, translated_output, output_csv)
    
    print(f"Reverse complement, mRNA sequence, and translation complete. Results saved to {reverse_complement_output}, {mrna_output}, and {translated_output}.")

# Usage
input_file = os.path.join(os.path.dirname(__file__), '../Sample_Data/coronavirus.fasta')
output_fasta = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/coronavirus_seq_manipulation.fasta')
output_csv = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/coronavirus_translation_stats.csv')
sequence_manipulation(input_file, output_fasta, output_csv)

