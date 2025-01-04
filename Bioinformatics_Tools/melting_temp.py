import os
from Bio import SeqIO
import csv
def check_empty_sequence(seq):
    if len(seq) == 0:
        print("Empty sequence detected.") #Debug print statement
        return True
    return False


def calculate_gc_content(seq):
    if check_empty_sequence(seq):  
        return 0.0
    g_count=seq.count('G')
    c_count=seq.count('C')
    return ((g_count+c_count)/len(seq))*100


def calculate_at_content(seq):
    if check_empty_sequence(seq):  
        return 0.0
    a_count=seq.count('A')
    t_count=seq.count('T')
    return ((a_count+t_count)/len(seq))*100

def calculating_melting_temp(seq):
    if check_empty_sequence(seq):  
        return 0.0
    gc_content=calculate_gc_content(seq)
    at_content=calculate_at_content(seq)
    return (2*(at_content)*len(seq)/100+4*(gc_content)*len(seq)/100)

def process_fasta_file(input_file,output_file=None):
    results=[]
    for record in SeqIO.parse(input_file, "fasta"):
        seq_id=record.id 
        seq=str(record.seq).upper()
        gc_content=calculate_gc_content(seq)
        at_content=calculate_at_content(seq)
        temp=calculating_melting_temp(seq)
        results.append([record.id,gc_content,at_content,temp])
        print(f"Sequence ID: {seq_id}")
        print(f"GC content: {gc_content:.2f}%")
        print(f"AT content: {at_content:.2f}%")
        print(f"Melting temperature: {temp:.2f}")
    if output_file:
        output_dir = os.path.join(os.path.dirname(input_file), 'output')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_file_path = os.path.join(output_dir, output_file)
        with open(output_file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Sequence ID', 'GC content', 'AT content', 'Melting temperature'])
            writer.writerows(results)
        print(f"Results saved to {output_file_path}")

#Analuzing coronavirus.fasta from Sample_Data folder
input_file = os.path.join(os.path.dirname(__file__), '../Sample_Data/coronavirus.fasta')
output_file = os.path.join(os.path.dirname(__file__), '../Sample_Data/tm.csv')
process_fasta_file(input_file, output_file)

    