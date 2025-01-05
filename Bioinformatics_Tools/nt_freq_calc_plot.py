'''                      Calculates nucleotide frequencies and plots them on a bar graph - works for any FASTA file with multiple sequences         '''
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import os
import numpy as np
import csv
from melting_temp import check_empty_sequence

# Calculate nucleotide frequencies
def calculate_freq(seq):
    if check_empty_sequence(seq):  
        return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}  # Return 0% if sequence is empty
    
    counter = Counter(seq)
    total_length = len(seq)

    return {
        'A': (counter['A'] / total_length) * 100,
        'T': (counter['T'] / total_length) * 100,
        'G': (counter['G'] / total_length) * 100,
        'C': (counter['C'] / total_length) * 100
    }

# Output nucleotide frequencies to CSV
def save_nt_freq_to_csv(frequencies, output_file):
    headers = ['Sequence ID', 'A Frequency (%)', 'T Frequency (%)', 'G Frequency (%)', 'C Frequency (%)']
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for seq_id, freq in frequencies.items():
            writer.writerow([seq_id] + list(freq.values()))
    print(f"Nucleotide frequencies saved to {output_file}")

# Visualize nucleotide frequencies as a grouped bar chart
def plot_nt_freq(frequencies, output_image, show_plot=False):
    plt.figure(figsize=(10, 6))
    nucleotides=['A', 'T', 'G', 'C']
    x_positions = np.arange(4) # A, T, G, C
    bar_width = 0.1  

    # Generate a color palette for different sequences in the FASTA file
    colors = plt.cm.tab20(np.linspace(0, 1, len(frequencies)))

    for idx, (seq_id, freq) in enumerate(frequencies.items()):
        freq_values = [freq[nt] for nt in nucleotides]  # Extract frequencies in A, T, G, C order
        offset = idx * bar_width - (bar_width * len(frequencies) / 2) 
        plt.bar(x_positions + offset, freq_values, bar_width, label=seq_id, color=colors[idx])

    plt.xlabel('Nucleotide', fontsize=12)
    plt.ylabel('Frequency (%)', fontsize=12)
    plt.title('Nucleotide Composition', fontsize=14)
    plt.xticks(x_positions, nucleotides, fontsize=10)
    plt.legend(title="Sequence ID", loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
    plt.tight_layout()
    
    plt.savefig(output_image)
    if show_plot:
        plt.show()
    plt.close()
    print(f"Grouped bar plot saved to {output_image}")

# Main function to process FASTA file, calculate frequencies, and generate outputs
def process_fasta_file(input_file, output_csv=None, output_image=None):
    nt_frequencies = {}

    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq).upper()
        
        # Calculate frequencies for each sequence
        nt_frequencies[seq_id] = calculate_freq(seq)

        print(f"Sequence ID: {seq_id}")
        for nt, freq in nt_frequencies[seq_id].items():
            print(f"{nt} frequency: {freq:.2f}%")
    
    if output_csv:
        save_nt_freq_to_csv(nt_frequencies, output_csv)
    if output_image:
        plot_nt_freq(nt_frequencies, output_image)
# Processing coronavirus.fasta
input_file = os.path.join(os.path.dirname(__file__), '../Sample_Data/coronavirus.fasta')
output_csv = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/nucleotide_frequencies.csv')
output_image = os.path.join(os.path.dirname(__file__), '../Sample_Data/output/nucleotide_frequencies_plot.png')
process_fasta_file(input_file, output_csv=output_csv, output_image=output_image)
