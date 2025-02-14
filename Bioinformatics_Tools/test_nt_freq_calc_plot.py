import pytest
import os
import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import write
from nt_freq_calc_plot import calculate_freq, save_nt_freq_to_csv, plot_nt_freq

# Test `calculate_freq` function
@pytest.mark.parametrize(
    "sequence, expected_frequencies",
    [
        ("AATTGGCC", {'A': 25.0, 'T': 25.0, 'G': 25.0, 'C': 25.0}),
        ("AAAA", {'A': 100.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}),
        ("", {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}),
    ],
)
def test_calculate_freq(sequence, expected_frequencies):
    assert calculate_freq(sequence) == expected_frequencies


# Test `save_nt_freq_to_csv` function
def test_save_nt_freq_to_csv(tmp_path):
    frequencies = {
        "seq1": {'A': 25.0, 'T': 25.0, 'G': 25.0, 'C': 25.0},
        "seq2": {'A': 50.0, 'T': 0.0, 'G': 25.0, 'C': 25.0},
    }
    output_file = tmp_path / "test_output.csv"

    save_nt_freq_to_csv(frequencies, output_file)

    with open(output_file, 'r') as f:
        reader = csv.reader(f)
        rows = list(reader)

    assert rows[0] == ['Sequence ID', 'A Frequency (%)', 'T Frequency (%)', 'G Frequency (%)', 'C Frequency (%)']
    assert rows[1] == ['seq1', '25.0', '25.0', '25.0', '25.0']
    assert rows[2] == ['seq2', '50.0', '0.0', '25.0', '25.0']


# Test `plot_nt_freq` function
def test_plot_nt_freq(tmp_path):
    frequencies = {
        "seq1": {'A': 25.0, 'T': 25.0, 'G': 25.0, 'C': 25.0},
        "seq2": {'A': 50.0, 'T': 0.0, 'G': 25.0, 'C': 25.0},
    }
    output_image = tmp_path / "test_plot.png"
    plot_nt_freq(frequencies, output_image)

    assert os.path.exists(output_image)


# Test complete pipeline with a FASTA file
def test_pipeline_with_fasta(tmp_path):
    # Create a sample FASTA file
    fasta_file = tmp_path / "test.fasta"
    records = [
        SeqRecord(Seq("AATTGGCC"), id="seq1", description=""),
        SeqRecord(Seq("AAAATTTT"), id="seq2", description=""),
    ]
    with open(fasta_file, 'w') as f:
        write(records, f, "fasta")

    # Test the full process
    from nt_freq_calc_plot import process_fasta_file
    output_csv = tmp_path / "test_output.csv"
    output_image = tmp_path / "test_plot.png"

    process_fasta_file(fasta_file, output_csv=output_csv, output_image=output_image)

    assert os.path.exists(output_csv)
    assert os.path.exists(output_image)

