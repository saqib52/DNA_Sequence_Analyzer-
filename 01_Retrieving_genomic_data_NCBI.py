from Bio import Entrez

# Set the email to be used with NCBI Entrez
Entrez.email = "saqibjawad52@gmail.com"  # Your email address

def fetch_fasta_sequence(accession_number):
    # Use the Entrez module to search and fetch the sequence
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    fasta_sequence = handle.read()
    handle.close()
    
    return fasta_sequence

# Example accession number for a genomic sequence (human chromosome 1)
accession_number = "NM_001301717.2"  # A smaller Gene sequence 
fasta_sequence = fetch_fasta_sequence(accession_number)

# Print or save the retrieved FASTA sequence
print(fasta_sequence)
