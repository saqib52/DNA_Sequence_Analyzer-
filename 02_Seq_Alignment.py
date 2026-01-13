# Using Biopython to run BLAST: You can use Biopython’s NCBIWWW module to submit your sequence to NCBI's BLAST server. Example

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Your sequence (use the sequence you fetched earlier)
sequence = """
>NM_001301717.2 Homo sapiens gene
ATGGAGGAGGCGGAGGAGGGAGGAGGAGAGGAGGAAGGAGGAGGAAGGAGGAGGAAGGAGGAGGAGGA...
"""

# Perform BLAST search against the nucleotide database
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

# Parse the BLAST results
blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    print("BLAST results for:", blast_record.query)
    for alignment in blast_record.alignments:
        print(f"Alignment title: {alignment.title}")
        for hsp in alignment.hsps:
            print(f"  Score: {hsp.score}, E-value: {hsp.expect}")
