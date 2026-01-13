import streamlit as st
import pandas as pd
import io
from collections import Counter
import re

# Page configuration
st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

# Title and description
st.title("🧬 DNA Sequence Analysis Tool")
st.markdown("Comprehensive DNA sequence analysis with primer design, restriction sites, motif search, visualizations, SNP detection, and export capabilities.")

# Sidebar for input
with st.sidebar:
    st.header("Input Options")
    input_method = st.radio(
        "Choose input method:",
        ["Text Input", "File Upload"]
    )
    
    if input_method == "Text Input":
        sequence_input = st.text_area(
            "Enter DNA Sequence:",
            height=150,
            placeholder="Enter DNA sequence here (e.g., ATGCGATCGATCG...)"
        )
    else:
        uploaded_file = st.file_uploader(
            "Upload FASTA file",
            type=['fasta', 'fa', 'txt']
        )
        if uploaded_file:
            file_content = uploaded_file.read().decode("utf-8")
            # Parse FASTA format
            lines = file_content.strip().split('\n')
            sequence_input = ''.join([line.strip() for line in lines if not line.startswith('>')])
        else:
            sequence_input = ""

# Functions
def validate_dna(sequence):
    """Validate if sequence contains only valid DNA nucleotides"""
    sequence = sequence.upper().replace(' ', '').replace('\n', '')
    valid_nucleotides = set('ATCGN')
    if not all(nuc in valid_nucleotides for nuc in sequence):
        return False, sequence
    return True, sequence

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence) - sequence.count('N')  # Exclude ambiguous bases
    if total == 0:
        return 0.0
    return (gc_count / total) * 100

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    sequence = sequence.upper()
    rev_comp = ''.join([complement.get(base, base) for base in reversed(sequence)])
    return rev_comp

def transcribe(sequence):
    """Transcribe DNA to RNA (replace T with U)"""
    sequence = sequence.upper()
    rna = sequence.replace('T', 'U')
    return rna

def translate(sequence):
    """Translate RNA sequence to protein using standard genetic code"""
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    sequence = sequence.upper()
    protein = ""
    # Find start codon (AUG)
    start_pos = sequence.find('AUG')
    if start_pos == -1:
        return "No start codon (AUG) found"
    
    # Translate from start codon
    for i in range(start_pos, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = genetic_code.get(codon, 'X')
            protein += amino_acid
            if amino_acid == '*':  # Stop codon
                break
    
    return protein

def get_sequence_statistics(sequence):
    """Calculate comprehensive sequence statistics"""
    sequence = sequence.upper()
    stats = {
        'length': len(sequence),
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C'),
        'N': sequence.count('N'),  # Ambiguous bases
    }
    stats['GC_content'] = calculate_gc_content(sequence)
    stats['AT_content'] = 100 - stats['GC_content']
    return stats

# Primer Design Functions
def calculate_tm(sequence, method='basic'):
    """Calculate melting temperature using simplified formula"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    length = len(sequence)
    
    if method == 'basic' and length <= 14:
        # For short sequences: Tm = 4(G+C) + 2(A+T)
        tm = 4 * gc_count + 2 * at_count
    else:
        # For longer sequences: Tm = 64.9 + 41*(G+C-16.4)/(A+T+G+C)
        if length > 0:
            tm = 64.9 + 41 * (gc_count - 16.4) / length
        else:
            tm = 0
    
    return round(tm, 2)

def check_primer_quality(primer):
    """Basic primer quality checks"""
    primer = primer.upper()
    issues = []
    
    # Check length
    if len(primer) < 18:
        issues.append("Too short (recommended: 18-25 bp)")
    elif len(primer) > 30:
        issues.append("Too long (recommended: 18-25 bp)")
    
    # Check GC content
    gc = calculate_gc_content(primer)
    if gc < 40:
        issues.append("Low GC content (< 40%)")
    elif gc > 60:
        issues.append("High GC content (> 60%)")
    
    # Check for runs
    if 'AAAA' in primer or 'TTTT' in primer or 'GGGG' in primer or 'CCCC' in primer:
        issues.append("Contains runs of 4+ identical bases")
    
    # Check 3' end (should end with G or C)
    if primer[-1] not in 'GC':
        issues.append("3' end should preferably be G or C")
    
    return issues

def design_primers(sequence, start_pos, end_pos, primer_length=20):
    """Design forward and reverse primers for a sequence region"""
    if end_pos > len(sequence) or start_pos < 0 or start_pos >= end_pos:
        return None, None
    
    # Forward primer (from start)
    forward = sequence[start_pos:start_pos + primer_length]
    
    # Reverse primer (reverse complement of end region)
    rev_start = max(0, end_pos - primer_length)
    rev_seq = sequence[rev_start:end_pos]
    reverse = reverse_complement(rev_seq)
    
    return forward, reverse

# Restriction Enzyme Database
RESTRICTION_ENZYMES = {
    'EcoRI': ('GAATTC', 'G^AATTC'),
    'BamHI': ('GGATCC', 'G^GATCC'),
    'HindIII': ('AAGCTT', 'A^AGCTT'),
    'XhoI': ('CTCGAG', 'C^TCGAG'),
    'NotI': ('GCGGCCGC', 'GC^GGCCGC'),
    'SalI': ('GTCGAC', 'G^TCGAC'),
    'PstI': ('CTGCAG', 'CTGCA^G'),
    'SmaI': ('CCCGGG', 'CCC^GGG'),
    'KpnI': ('GGTACC', 'GGTAC^C'),
    'SacI': ('GAGCTC', 'GAGCT^C'),
    'NdeI': ('CATATG', 'CA^TATG'),
    'XbaI': ('TCTAGA', 'T^CTAGA'),
    'SpeI': ('ACTAGT', 'A^CTAGT'),
    'SphI': ('GCATGC', 'GCATG^C'),
    'NcoI': ('CCATGG', 'C^CATGG'),
}

def find_restriction_sites(sequence):
    """Find restriction enzyme recognition sites in sequence"""
    sites = []
    sequence = sequence.upper()
    
    for enzyme, (recognition_site, _) in RESTRICTION_ENZYMES.items():
        pattern = recognition_site
        positions = []
        start = 0
        while True:
            pos = sequence.find(pattern, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        if positions:
            sites.append({
                'enzyme': enzyme,
                'recognition_site': recognition_site,
                'positions': positions,
                'count': len(positions)
            })
    
    return sorted(sites, key=lambda x: x['count'], reverse=True)

# Motif Search
def find_motif(sequence, motif, case_sensitive=False):
    """Find all occurrences of a motif in the sequence"""
    if not case_sensitive:
        sequence = sequence.upper()
        motif = motif.upper()
    
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    
    return positions

def find_motif_regex(sequence, pattern, case_sensitive=False):
    """Find motifs using regex pattern"""
    if not case_sensitive:
        flags = re.IGNORECASE
    else:
        flags = 0
    
    matches = []
    for match in re.finditer(pattern, sequence, flags):
        matches.append({
            'start': match.start(),
            'end': match.end(),
            'sequence': match.group()
        })
    
    return matches

# Codon Frequency Analysis
def get_codon_frequency(sequence):
    """Calculate codon frequency from DNA sequence"""
    sequence = sequence.upper().replace('T', 'U')  # Convert to RNA
    codons = []
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    
    return Counter(codons)

# SNP Detection
def detect_snps(seq1, seq2):
    """Detect single nucleotide polymorphisms between two sequences"""
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    snps = []
    min_len = min(len(seq1), len(seq2))
    
    for i in range(min_len):
        if seq1[i] != seq2[i] and seq1[i] in 'ATCG' and seq2[i] in 'ATCG':
            snps.append({
                'position': i + 1,  # 1-indexed
                'ref': seq1[i],
                'alt': seq2[i],
                'type': 'SNP'
            })
    
    return snps

# GC Content Plot Data
def gc_content_plot_data(sequence, window_size=50):
    """Calculate GC content in sliding windows"""
    sequence = sequence.upper()
    gc_values = []
    positions = []
    
    for i in range(0, len(sequence) - window_size + 1, window_size // 2):
        window = sequence[i:i + window_size]
        gc = calculate_gc_content(window)
        gc_values.append(gc)
        positions.append(i + window_size // 2)
    
    return positions, gc_values

# Export Functions
def export_to_fasta(sequence, name="sequence"):
    """Export sequence to FASTA format"""
    return f">{name}\n{sequence}\n"

def export_to_csv(data, filename="export"):
    """Export data to CSV format"""
    df = pd.DataFrame(data)
    csv = df.to_csv(index=False)
    return csv

def export_results_text(results_dict):
    """Export results as formatted text"""
    output = []
    for key, value in results_dict.items():
        output.append(f"{key}: {value}\n")
    return "".join(output)

# Main analysis
if sequence_input:
    # Validate and clean sequence
    is_valid, cleaned_sequence = validate_dna(sequence_input)
    
    if not cleaned_sequence:
        st.warning("Please enter a DNA sequence.")
    elif not is_valid:
        st.error("❌ Invalid sequence! Please enter only A, T, G, C, or N (for ambiguous bases).")
    else:
        st.success(f"✅ Valid DNA sequence loaded! Length: {len(cleaned_sequence)} bp")
        
        # Create tabs for different analyses
        tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10 = st.tabs([
            "📊 Statistics", 
            "🧪 GC Content", 
            "🔄 Reverse Complement", 
            "📝 Transcription", 
            "🔄 Translation",
            "🧬 Primer Design",
            "✂️ Restriction Sites",
            "🔍 Motif Search",
            "📈 Visualizations",
            "🔬 SNP Detection"
        ])
        
        with tab1:
            st.header("Sequence Statistics")
            stats = get_sequence_statistics(cleaned_sequence)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Length", f"{stats['length']} bp")
                st.metric("GC Content", f"{stats['GC_content']:.2f}%")
                st.metric("AT Content", f"{stats['AT_content']:.2f}%")
            
            with col2:
                st.metric("Adenine (A)", stats['A'])
                st.metric("Thymine (T)", stats['T'])
                st.metric("Guanine (G)", stats['G'])
            
            with col3:
                st.metric("Cytosine (C)", stats['C'])
                if stats['N'] > 0:
                    st.metric("Ambiguous (N)", stats['N'])
            
            # Nucleotide composition chart
            if stats['N'] == 0:
                st.subheader("Nucleotide Composition")
                chart_data = {
                    'A': stats['A'],
                    'T': stats['T'],
                    'G': stats['G'],
                    'C': stats['C']
                }
                st.bar_chart(chart_data)
        
        with tab2:
            st.header("GC Content Analysis")
            gc_content = calculate_gc_content(cleaned_sequence)
            
            col1, col2 = st.columns([1, 2])
            with col1:
                st.metric("GC Content", f"{gc_content:.2f}%")
                
                # GC content interpretation
                if gc_content < 35:
                    st.info("Low GC content (< 35%)")
                elif gc_content > 65:
                    st.info("High GC content (> 65%)")
                else:
                    st.info("Moderate GC content (35-65%)")
            
            with col2:
                st.subheader("GC Content Visualization")
                gc_data = {
                    'GC': gc_content,
                    'AT': 100 - gc_content
                }
                st.bar_chart(gc_data)
        
        with tab3:
            st.header("Reverse Complement")
            rev_comp = reverse_complement(cleaned_sequence)
            
            st.subheader("Original Sequence")
            st.code(cleaned_sequence, language=None)
            
            st.subheader("Reverse Complement")
            st.code(rev_comp, language=None)
            
            # Copy button functionality
            st.text_area(
                "Copy reverse complement:",
                rev_comp,
                height=100,
                key="rev_comp_output",
                disabled=False
            )
        
        with tab4:
            st.header("Transcription (DNA → RNA)")
            rna_sequence = transcribe(cleaned_sequence)
            
            st.subheader("DNA Sequence")
            st.code(cleaned_sequence, language=None)
            
            st.subheader("RNA Sequence (Transcribed)")
            st.code(rna_sequence, language=None)
            
            st.info("💡 Transcription converts T (Thymine) to U (Uracil)")
            
            st.text_area(
                "Copy RNA sequence:",
                rna_sequence,
                height=100,
                key="rna_output",
                disabled=False
            )
        
        with tab5:
            st.header("Translation (RNA → Protein)")
            
            # First transcribe to RNA
            rna_seq = transcribe(cleaned_sequence)
            
            st.subheader("RNA Sequence")
            st.code(rna_seq, language=None)
            
            st.subheader("Protein Sequence (Translated)")
            protein_seq = translate(rna_seq)
            
            if protein_seq.startswith("No start"):
                st.warning(protein_seq)
            else:
                st.code(protein_seq, language=None)
                st.metric("Protein Length", f"{len(protein_seq)} amino acids")
                
                st.text_area(
                    "Copy protein sequence:",
                    protein_seq,
                    height=100,
                    key="protein_output",
                    disabled=False
                )
                
                st.info("💡 Translation uses standard genetic code. Start codon: AUG (M), Stop codons: UAA, UAG, UGA (*)")
        
        with tab6:
            st.header("🧬 Primer Design")
            
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Manual Primer Entry")
                forward_primer = st.text_input("Forward Primer:", placeholder="Enter forward primer sequence")
                reverse_primer = st.text_input("Reverse Primer:", placeholder="Enter reverse primer sequence")
                
                if forward_primer or reverse_primer:
                    if forward_primer:
                        st.write("**Forward Primer Analysis:**")
                        f_gc = calculate_gc_content(forward_primer)
                        f_tm = calculate_tm(forward_primer)
                        f_issues = check_primer_quality(forward_primer)
                        
                        col_f1, col_f2 = st.columns(2)
                        with col_f1:
                            st.metric("GC Content", f"{f_gc:.2f}%")
                            st.metric("Tm", f"{f_tm}°C")
                        with col_f2:
                            st.write("**Sequence:**")
                            st.code(forward_primer, language=None)
                        
                        if f_issues:
                            st.warning("⚠️ Issues found:")
                            for issue in f_issues:
                                st.write(f"- {issue}")
                        else:
                            st.success("✅ Primer quality is good!")
                    
                    if reverse_primer:
                        st.write("**Reverse Primer Analysis:**")
                        r_gc = calculate_gc_content(reverse_primer)
                        r_tm = calculate_tm(reverse_primer)
                        r_issues = check_primer_quality(reverse_primer)
                        
                        col_r1, col_r2 = st.columns(2)
                        with col_r1:
                            st.metric("GC Content", f"{r_gc:.2f}%")
                            st.metric("Tm", f"{r_tm}°C")
                        with col_r2:
                            st.write("**Sequence:**")
                            st.code(reverse_primer, language=None)
                        
                        if r_issues:
                            st.warning("⚠️ Issues found:")
                            for issue in r_issues:
                                st.write(f"- {issue}")
                        else:
                            st.success("✅ Primer quality is good!")
                    
                    if forward_primer and reverse_primer:
                        tm_diff = abs(f_tm - r_tm)
                        if tm_diff > 5:
                            st.warning(f"⚠️ Large Tm difference: {tm_diff:.2f}°C (recommended: < 5°C)")
                        else:
                            st.success(f"✅ Good Tm match: {tm_diff:.2f}°C difference")
            
            with col2:
                st.subheader("Automatic Primer Design")
                st.write(f"Sequence length: {len(cleaned_sequence)} bp")
                start_pos = st.number_input("Start position:", min_value=0, max_value=len(cleaned_sequence)-1, value=0)
                end_pos = st.number_input("End position:", min_value=1, max_value=len(cleaned_sequence), value=min(200, len(cleaned_sequence)))
                primer_len = st.slider("Primer length:", min_value=15, max_value=30, value=20)
                
                if st.button("Design Primers"):
                    if end_pos > start_pos:
                        forward, reverse = design_primers(cleaned_sequence, start_pos, end_pos, primer_len)
                        if forward and reverse:
                            st.success("Primers designed successfully!")
                            
                            st.write("**Designed Primers:**")
                            st.write(f"**Forward:** {forward}")
                            st.write(f"**Reverse:** {reverse}")
                            
                            f_tm = calculate_tm(forward)
                            r_tm = calculate_tm(reverse)
                            st.write(f"Forward Tm: {f_tm}°C | Reverse Tm: {r_tm}°C")
                        else:
                            st.error("Invalid primer design parameters")
        
        with tab7:
            st.header("✂️ Restriction Enzyme Sites")
            
            sites = find_restriction_sites(cleaned_sequence)
            
            if sites:
                st.write(f"Found {len(sites)} enzyme(s) with recognition sites:")
                
                # Summary table
                summary_data = []
                for site in sites:
                    summary_data.append({
                        'Enzyme': site['enzyme'],
                        'Recognition Site': site['recognition_site'],
                        'Count': site['count'],
                        'Positions': ', '.join([str(p+1) for p in site['positions'][:5]]) + ('...' if len(site['positions']) > 5 else '')
                    })
                
                df_sites = pd.DataFrame(summary_data)
                st.dataframe(df_sites, use_container_width=True)
                
                # Detailed view
                selected_enzyme = st.selectbox("Select enzyme for detailed view:", [s['enzyme'] for s in sites])
                if selected_enzyme:
                    enzyme_info = next(s for s in sites if s['enzyme'] == selected_enzyme)
                    st.write(f"**{selected_enzyme}** - Recognition site: {enzyme_info['recognition_site']}")
                    st.write(f"**Number of sites:** {enzyme_info['count']}")
                    st.write(f"**Positions (1-indexed):** {', '.join([str(p+1) for p in enzyme_info['positions']])}")
            else:
                st.info("No restriction enzyme recognition sites found in the sequence.")
            
            st.info("💡 Database includes: EcoRI, BamHI, HindIII, XhoI, NotI, SalI, PstI, SmaI, KpnI, SacI, NdeI, XbaI, SpeI, SphI, NcoI")
        
        with tab8:
            st.header("🔍 Motif Search")
            
            search_type = st.radio("Search type:", ["Exact match", "Regex pattern"], horizontal=True)
            
            if search_type == "Exact match":
                motif = st.text_input("Enter motif to search:", placeholder="e.g., ATGC")
                case_sensitive = st.checkbox("Case sensitive", value=False)
                
                if motif:
                    positions = find_motif(cleaned_sequence, motif, case_sensitive)
                    if positions:
                        st.success(f"Found {len(positions)} occurrence(s)")
                        st.write("**Positions (1-indexed):**", ", ".join([str(p+1) for p in positions]))
                        
                        # Show sequences around matches
                        st.subheader("Matches:")
                        for i, pos in enumerate(positions[:10]):  # Show first 10
                            start = max(0, pos - 10)
                            end = min(len(cleaned_sequence), pos + len(motif) + 10)
                            context = cleaned_sequence[start:end]
                            highlighted = context.replace(motif.upper(), f"**{motif.upper()}**")
                            st.write(f"{pos+1}: ...{context}...")
                    else:
                        st.warning("Motif not found in sequence")
            else:
                pattern = st.text_input("Enter regex pattern:", placeholder="e.g., [AT]{3,5}")
                case_sensitive = st.checkbox("Case sensitive", value=False)
                
                if pattern:
                    try:
                        matches = find_motif_regex(cleaned_sequence, pattern, case_sensitive)
                        if matches:
                            st.success(f"Found {len(matches)} match(es)")
                            
                            match_data = []
                            for match in matches[:20]:  # Show first 20
                                match_data.append({
                                    'Start': match['start'] + 1,
                                    'End': match['end'],
                                    'Length': match['end'] - match['start'],
                                    'Sequence': match['sequence']
                                })
                            
                            df_matches = pd.DataFrame(match_data)
                            st.dataframe(df_matches, use_container_width=True)
                        else:
                            st.warning("No matches found")
                    except re.error as e:
                        st.error(f"Invalid regex pattern: {e}")
        
        with tab9:
            st.header("📈 Interactive Visualizations")
            
            viz_type = st.selectbox("Select visualization:", 
                                   ["GC Content Plot", "Nucleotide Frequency", "Codon Frequency"])
            
            if viz_type == "GC Content Plot":
                window_size = st.slider("Window size:", min_value=20, max_value=200, value=50, step=10)
                positions, gc_values = gc_content_plot_data(cleaned_sequence, window_size)
                
                if positions:
                    gc_df = pd.DataFrame({
                        'Position': positions,
                        'GC Content (%)': gc_values
                    })
                    st.line_chart(gc_df.set_index('Position'))
                    st.caption(f"GC content in sliding windows of {window_size} bp")
            
            elif viz_type == "Nucleotide Frequency":
                stats = get_sequence_statistics(cleaned_sequence)
                freq_data = pd.DataFrame({
                    'Nucleotide': ['A', 'T', 'G', 'C'],
                    'Count': [stats['A'], stats['T'], stats['G'], stats['C']],
                    'Percentage': [
                        (stats['A']/stats['length']*100) if stats['length'] > 0 else 0,
                        (stats['T']/stats['length']*100) if stats['length'] > 0 else 0,
                        (stats['G']/stats['length']*100) if stats['length'] > 0 else 0,
                        (stats['C']/stats['length']*100) if stats['length'] > 0 else 0
                    ]
                })
                st.bar_chart(freq_data.set_index('Nucleotide')['Count'])
                st.dataframe(freq_data, use_container_width=True)
            
            elif viz_type == "Codon Frequency":
                codon_freq = get_codon_frequency(cleaned_sequence)
                
                if codon_freq:
                    codon_df = pd.DataFrame({
                        'Codon': list(codon_freq.keys()),
                        'Frequency': list(codon_freq.values())
                    }).sort_values('Frequency', ascending=False)
                    
                    st.bar_chart(codon_df.set_index('Codon')['Frequency'])
                    st.dataframe(codon_df.head(20), use_container_width=True)
                    st.caption("Top 20 most frequent codons")
        
        with tab10:
            st.header("🔬 SNP Detection")
            
            st.write("Compare two sequences to detect Single Nucleotide Polymorphisms (SNPs)")
            
            comparison_method = st.radio("Input method:", ["Paste second sequence", "Upload FASTA file"], horizontal=True)
            
            seq2_input = ""
            if comparison_method == "Paste second sequence":
                seq2_input = st.text_area("Enter second sequence for comparison:", height=150)
            else:
                uploaded_file2 = st.file_uploader("Upload second sequence (FASTA)", type=['fasta', 'fa', 'txt'], key="snp_file")
                if uploaded_file2:
                    file_content2 = uploaded_file2.read().decode("utf-8")
                    lines2 = file_content2.strip().split('\n')
                    seq2_input = ''.join([line.strip() for line in lines2 if not line.startswith('>')])
            
            if seq2_input:
                is_valid2, cleaned_seq2 = validate_dna(seq2_input)
                
                if is_valid2 and cleaned_seq2:
                    snps = detect_snps(cleaned_sequence, cleaned_seq2)
                    
                    if snps:
                        st.success(f"Found {len(snps)} SNP(s)")
                        
                        snp_data = []
                        for snp in snps:
                            snp_data.append({
                                'Position': snp['position'],
                                'Reference': snp['ref'],
                                'Alternate': snp['alt'],
                                'Change': f"{snp['ref']} → {snp['alt']}"
                            })
                        
                        df_snps = pd.DataFrame(snp_data)
                        st.dataframe(df_snps, use_container_width=True)
                        
                        # Summary
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Total SNPs", len(snps))
                        with col2:
                            transitions = sum(1 for s in snps if (s['ref'] in 'AG' and s['alt'] in 'AG') or (s['ref'] in 'CT' and s['alt'] in 'CT'))
                            st.metric("Transitions", transitions)
                        with col3:
                            transversions = len(snps) - transitions
                            st.metric("Transversions", transversions)
                    else:
                        st.info("No SNPs detected - sequences are identical")
                else:
                    st.error("Invalid second sequence")
        
        # Export section
        st.divider()
        with st.expander("💾 Export Results"):
            export_format = st.selectbox("Export format:", ["FASTA", "CSV (Statistics)", "TXT (Summary)"])
            
            if export_format == "FASTA":
                seq_name = st.text_input("Sequence name:", value="sequence")
                export_data = export_to_fasta(cleaned_sequence, seq_name)
                st.download_button(
                    label="Download FASTA",
                    data=export_data,
                    file_name=f"{seq_name}.fasta",
                    mime="text/plain"
                )
            
            elif export_format == "CSV (Statistics)":
                stats = get_sequence_statistics(cleaned_sequence)
                export_data = export_to_csv([stats])
                st.download_button(
                    label="Download CSV",
                    data=export_data,
                    file_name="sequence_statistics.csv",
                    mime="text/csv"
                )
            
            else:
                stats = get_sequence_statistics(cleaned_sequence)
                export_data = export_results_text(stats)
                st.download_button(
                    label="Download TXT",
                    data=export_data,
                    file_name="sequence_summary.txt",
                    mime="text/plain"
                )
        
        # Additional features in expander
        with st.expander("📋 Sequence Information"):
            st.write(f"**Sequence Length:** {len(cleaned_sequence)} nucleotides")
            st.write(f"**First 50 bases:** {cleaned_sequence[:50]}...")
            st.write(f"**Last 50 bases:** ...{cleaned_sequence[-50:]}")
            
else:
    st.info("👈 Please enter a DNA sequence in the sidebar to begin analysis.")
    
    # Example sequence
    with st.expander("📝 Example DNA Sequence"):
        example_seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        st.code(example_seq, language=None)
        st.caption("Copy and paste this sequence into the sidebar to try the tool!")
