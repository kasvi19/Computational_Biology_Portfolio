from sequence_manipulation import reverse_complement, mrna, translation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os
import tempfile
import pytest

# Define sample sequences for testing
@pytest.fixture
def sample_records():
    records = [
        SeqRecord(Seq("ATGCGTACG"), id="seq1"),
        SeqRecord(Seq("GGCTAATGC"), id="seq2"),
        SeqRecord(Seq("TACGGT"), id="seq3"),
    ]
    return records

# Test reverse complement function
def test_reverse_complement(sample_records):
    expected_rc_seqs = ["CGTACGCAT", "GCATTAGCC", "ACCGTA"]
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_output:
        reverse_complement(sample_records, temp_output.name)
        reverse_complemented_records = list(SeqIO.parse(temp_output.name, "fasta"))
    for record, expected_seq in zip(reverse_complemented_records, expected_rc_seqs):
        assert str(record.seq) == expected_seq

# Test mRNA function
def test_mrna(sample_records):
    expected_mrna_seqs = ["AUGCGUACG", "GGCUAAUGC", "UACGGU"]
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_output:
        mrna(sample_records, temp_output.name, is_template_strand=False)
        mrna_records = list(SeqIO.parse(temp_output.name, "fasta"))
    for record, expected_seq in zip(mrna_records, expected_mrna_seqs):
        assert str(record.seq) == expected_seq

# Test translation function
def test_translation(sample_records):
    expected_translated_seqs = ["MRT", "G*C", "YG"]
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_input_fasta, \
        tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_output_fasta, \
        tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as temp_output_csv:

        # Write the sample records to the temp input FASTA
        SeqIO.write(sample_records, temp_input_fasta.name, "fasta")

        # Call the translation function
        translated_records, stats = translation(temp_input_fasta.name, temp_output_fasta.name, temp_output_csv.name)

        # Validate translated sequences
        for record, expected_seq in zip(translated_records, expected_translated_seqs):
            assert str(record.seq) == expected_seq

        # Validate stats
        assert stats["total_sequences"] == 3
        assert stats["total_stop_codons"] == 1
        assert stats["total_aa"] == 8
        assert stats["invalid_codons"] == 0

# Clean up temp files after tests
@pytest.fixture(scope="function", autouse=True)
def cleanup_temp_files():
    temp_files = []

    def register(file_path):
        temp_files.append(file_path)
        return file_path

    yield register

    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)