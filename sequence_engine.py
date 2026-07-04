import os
import re
import json
from Bio.Seq import Seq
import ml_filter

# Load reference database of AMR genes
REF_DB_PATH = os.path.join(os.path.dirname(__file__), "data", "amr_reference.json")

def load_references() -> list:
    """
    Loads reference AMR sequences from data/amr_reference.json.
    """
    if not os.path.exists(REF_DB_PATH):
        return []
    try:
        with open(REF_DB_PATH, "r") as f:
            return json.load(f)
    except Exception:
        return []

def sanitize_sequence(sequence: str) -> str:
    """
    Sanitizes raw sequence input. If strings contain elements outside 
    valid IUPAC nucleotide codes (A, T, G, C, N), execution must abort.
    Also strips FASTA headers and whitespace.
    """
    # Remove FASTA header if present
    cleaned = re.sub(r'^>.*?\n', '', sequence)
    # Remove any headers in the middle (if multiple sequences)
    cleaned = re.sub(r'>.*?\n', '', cleaned)
    # Remove all whitespace, newlines, and carriage returns
    cleaned = re.sub(r'\s+', '', cleaned).upper()
    
    if not cleaned:
        raise ValueError("Invalid Input Error: Sequence input is empty or invalid.")
        
    # Check for invalid IUPAC characters
    invalid_chars = set(cleaned) - {'A', 'T', 'G', 'C', 'N'}
    if invalid_chars:
        # Abort execution due to tampering/injection threat (STRIDE)
        raise ValueError(
            f"Security Sanitization Alert: Sequence contains invalid non-IUPAC nucleotide characters: "
            f"{sorted(list(invalid_chars))}. Only A, T, G, C, N are allowed."
        )
    return cleaned

def align_needleman_wunsch(seq1: str, seq2: str) -> float:
    """
    Computes global alignment between two sequences using Needleman-Wunsch algorithm.
    Returns normalized sequence identity (0.0 to 1.0).
    """
    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return 0.0
        
    match_score = 2
    mismatch_score = -1
    gap_penalty = -2
    
    # Initialize DP matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty
        
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match, delete, insert)
            
    # Traceback to calculate sequence identity
    i, j = m, n
    matches = 0
    align_len = 0
    while i > 0 or j > 0:
        align_len += 1
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            if seq1[i-1] == seq2[j-1]:
                matches += 1
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap_penalty:
            i -= 1
        else:
            j -= 1
            
    if align_len == 0:
        return 0.0
    return matches / min(m, n)

def extract_and_analyze_orfs(dna_sequence: str) -> list:
    """
    Translates DNA in all 6 reading frames, extracts ORFs >= 30 aa,
    runs the ML filter, checks against reference databases, and tags pseudogenes.
    """
    sanitized_dna = sanitize_sequence(dna_sequence)
    references = load_references()
    seq_obj = Seq(sanitized_dna)
    
    # Prepare the 6 frames
    frames = []
    # Forward frames (0, 1, 2)
    for frame in range(3):
        frames.append((seq_obj[frame:], frame, False))
    # Reverse frames (0, 1, 2)
    rev_seq = seq_obj.reverse_complement()
    for frame in range(3):
        frames.append((rev_seq[frame:], frame, True))
        
    extracted_orfs = []
    
    for s, frame, is_reverse in frames:
        # Translate the full reading frame, keeping stop codons as '*'
        translated = str(s.translate(to_stop=False))
        n = len(translated)
        
        # 1. Standard ORF Extraction: Scan codon by codon for 'M'
        i = 0
        while i < n:
            if translated[i] == 'M':
                # Find the next stop codon '*'
                stop_idx = translated.find('*', i)
                if stop_idx != -1:
                    aa_seq = translated[i:stop_idx]
                    next_i = stop_idx + 1
                else:
                    aa_seq = translated[i:]
                    next_i = n
                    
                # Biological noise filter: discard sequences under 30 aa
                if len(aa_seq) >= 30:
                    # Run ML filter to remove out-of-frame translation noise
                    is_coding = ml_filter.predict_is_protein(aa_seq)
                    if is_coding:
                        # Compare against references
                        best_match = None
                        best_identity = 0.0
                        
                        for ref in references:
                            identity = align_needleman_wunsch(aa_seq, ref["amino_acid_sequence"])
                            if identity > best_identity:
                                best_identity = identity
                                best_match = ref
                                
                        # Tag classification
                        status = "NOVEL_PROTEIN"
                        action_required = True
                        ref_name = None
                        
                        if best_match and best_identity >= 0.8:
                            status = "MUTATED_VARIANT"
                            ref_name = best_match["gene_name"]
                        
                        extracted_orfs.append({
                            "amino_acid_sequence": aa_seq,
                            "status": status,
                            "action_required": action_required,
                            "sequence_identity": round(best_identity, 3),
                            "reference_match": ref_name,
                            "frame": frame,
                            "is_reverse": is_reverse,
                            "is_pseudogene": False
                        })
                i = next_i
            else:
                i += 1
                
        # 2. Pseudogene Heuristic Loop:
        # Check if translating through stops yields a sequence that aligns well to a reference,
        # but contains premature stop codons that segment the protein.
        # We do this by scanning translations with internal stop codons.
        i = 0
        while i < n:
            if translated[i] == 'M':
                # Look ahead for a potential sequence that spans stop codons
                # We can trace up to 300 aa or the end of the frame
                limit = min(i + 350, n)
                potential_full_translation = translated[i:limit]
                
                # Check alignment of this full segment (replacing '*' with 'X' to align through stops)
                clean_full_translation = potential_full_translation.replace('*', 'X')
                
                best_ref = None
                best_full_identity = 0.0
                for ref in references:
                    full_identity = align_needleman_wunsch(clean_full_translation, ref["amino_acid_sequence"])
                    if full_identity > best_full_identity:
                        best_full_identity = full_identity
                        best_ref = ref
                
                # If there's a strong alignment (>= 70%) to a reference protein,
                # AND there are internal stop codons within the first 85% of the reference length,
                # we tag this as a Pseudogene.
                if best_ref and best_full_identity >= 0.7:
                    ref_seq = best_ref["amino_acid_sequence"]
                    # Find if there is an in-frame stop codon before the expected end
                    first_stop = potential_full_translation.find('*')
                    if first_stop != -1 and first_stop < int(len(ref_seq) * 0.85):
                        # Extract the truncated sequence up to the first stop codon
                        truncated_seq = potential_full_translation[:first_stop]
                        # Discard very short fragments
                        if len(truncated_seq) >= 15:
                            # Verify if we haven't already added this exact sequence
                            already_added = any(o["amino_acid_sequence"] == truncated_seq for o in extracted_orfs)
                            if not already_added:
                                extracted_orfs.append({
                                    "amino_acid_sequence": truncated_seq,
                                    "status": "MUTATED_VARIANT", # Complies with output schema targets
                                    "action_required": False,
                                    "sequence_identity": round(best_full_identity, 3),
                                    "reference_match": best_ref["gene_name"],
                                    "frame": frame,
                                    "is_reverse": is_reverse,
                                    "is_pseudogene": True,
                                    "notes": "Pseudogene detected: premature stop codon identified."
                                })
                
                # Move to next stop codon to avoid infinite loops or overlaps
                stop_idx = translated.find('*', i)
                if stop_idx != -1:
                    i = stop_idx + 1
                else:
                    break
            else:
                i += 1
                
    return extracted_orfs
