import os
import shutil
import logging
import unittest
import warnings

# Suppress all third-party deprecation and packaging warnings
warnings.filterwarnings("ignore")

# Silence app and httpx loggers during test execution
logging.getLogger("AMR_Pipeline").setLevel(logging.CRITICAL)
logging.getLogger("httpx").setLevel(logging.CRITICAL)
logging.getLogger("uvicorn").setLevel(logging.CRITICAL)

from fastapi.testclient import TestClient

import ml_filter
import sequence_engine
import structural_docking
from app import app

# Simple reverse translation table for test sequence generation
CODON_MAP = {
    'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGC',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAT',
    '*': 'TAA'
}

def aa_to_dna(aa_seq: str) -> str:
    return "".join(CODON_MAP.get(aa, 'ATG') for aa in aa_seq)

class TestAMRPipeline(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Ensure ML model is trained and cached
        cls.clf = ml_filter.get_classifier()
        
        # Clean existing affinity matrix if present to start clean
        if os.path.exists(structural_docking.AFFINITY_MATRIX_PATH):
            os.remove(structural_docking.AFFINITY_MATRIX_PATH)
            
    def test_01_ml_noise_filter(self):
        print("\n--- Testing ML Noise Filter ---")
        # Test positive sequence (typical protein sequence)
        pos_seq = "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILES"
        self.assertTrue(ml_filter.predict_is_protein(pos_seq))
        
        # Test negative sequence (short or garbage sequence)
        neg_seq = "ACACAC"
        self.assertFalse(ml_filter.predict_is_protein(neg_seq))
        
    def test_02_iupac_sanitization(self):
        print("\n--- Testing IUPAC Sanitization ---")
        # Valid sequence should pass and be cleaned
        valid_dna = ">TestSeq\nATGGCCGCTAAATGA\n"
        cleaned = sequence_engine.sanitize_sequence(valid_dna)
        self.assertEqual(cleaned, "ATGGCCGCTAAATGA")
        
        # Invalid sequence (with amino acid characters like 'L', 'P') should raise ValueError
        invalid_dna = "ATGGKKVLMVLSLLVILPAVASAT"
        with self.assertRaises(ValueError) as context:
            sequence_engine.sanitize_sequence(invalid_dna)
        self.assertIn("Security Sanitization Alert", str(context.exception))
        
    def test_03_orf_extraction_and_alignment(self):
        print("\n--- Testing ORF Extraction and Alignment ---")
        # Construct DNA sequence for a portion of blaTEM-1 (first 50 amino acids)
        ref_db = sequence_engine.load_references()
        bla_tem = ref_db[0]["amino_acid_sequence"]
        
        # We need an ORF >= 30 aa starting with M (ATG)
        aa_fragment = bla_tem[:60] # MSIQH...
        dna_fragment = aa_to_dna(aa_fragment) + "TAA" # add stop codon
        
        # Analyze the DNA fragment
        orfs = sequence_engine.extract_and_analyze_orfs(dna_fragment)
        
        self.assertGreaterEqual(len(orfs), 1)
        # Check that the extracted sequence is identical
        matching_orf = [o for o in orfs if o["amino_acid_sequence"] == aa_fragment]
        self.assertEqual(len(matching_orf), 1)
        self.assertEqual(matching_orf[0]["status"], "MUTATED_VARIANT")
        self.assertEqual(matching_orf[0]["reference_match"], "blaTEM-1")
        self.assertTrue(matching_orf[0]["action_required"])
        
    def test_04_pseudogene_detection(self):
        print("\n--- Testing Pseudogene Heuristic ---")
        # Take a reference gene and introduce a premature stop codon
        ref_db = sequence_engine.load_references()
        bla_tem = ref_db[0]["amino_acid_sequence"]
        
        # MSIQHFRVAL... (introduce * at position 20)
        mutated_aa = bla_tem[:20] + "*" + bla_tem[21:100]
        mutated_dna = aa_to_dna(mutated_aa) + "TAA"
        
        orfs = sequence_engine.extract_and_analyze_orfs(mutated_dna)
        
        # Check if we successfully tagged a pseudogene
        pseudogenes = [o for o in orfs if o["is_pseudogene"]]
        self.assertGreaterEqual(len(pseudogenes), 1)
        self.assertFalse(pseudogenes[0]["action_required"])
        self.assertEqual(pseudogenes[0]["reference_match"], "blaTEM-1")
        
    def test_05_docking_and_matrix_append(self):
        print("\n--- Testing Structural Docking Matrix ---")
        test_aa = "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILES"
        
        # Run compound screen
        results = structural_docking.screen_compounds(test_aa, ref_match="blaTEM-1", seq_identity=0.95)
        
        self.assertEqual(len(results), len(structural_docking.ANTIBIOTICS))
        self.assertTrue(os.path.exists(structural_docking.AFFINITY_MATRIX_PATH))
        
        # Check the highest affinity compound (should be Penicillin G or Ampicillin for blaTEM-1)
        best_comp = min(results, key=lambda x: x["binding_affinity_kcal_mol"])
        self.assertIn(best_comp["compound_name"], ["Penicillin G", "Ampicillin"])
        self.assertLess(best_comp["binding_affinity_kcal_mol"], -7.0)
        
    def test_06_fastapi_endpoints(self):
        print("\n--- Testing FastAPI Endpoints ---")
        client = TestClient(app)
        
        # 1. Test POST with valid sequence
        ref_db = sequence_engine.load_references()
        valid_dna = aa_to_dna(ref_db[0]["amino_acid_sequence"][:80]) + "TAA"
        
        response = client.post("/api/v1/analyze_genome", json={"fasta_data": valid_dna})
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertTrue(data["success"])
        self.assertGreaterEqual(len(data["results"]), 1)
        self.assertEqual(data["results"][0]["reference_match"], "blaTEM-1")
        self.assertGreater(len(data["results"][0]["docking_results"]), 0)
        
        # 2. Test POST with invalid sequence (Tampering injection test)
        invalid_dna = "ATGGKKVLMVLSLLVILPAVASATGTAAATMKVSTLLLLLCVSLALVAAATGCCCCGGGG"
        response = client.post("/api/v1/analyze_genome", json={"fasta_data": invalid_dna})
        self.assertEqual(response.status_code, 400)
        self.assertIn("Security Sanitization Alert", response.json()["detail"])
        
        # 3. Test POST with empty data
        response = client.post("/api/v1/analyze_genome", json={"fasta_data": ""})
        self.assertEqual(response.status_code, 400)
        self.assertIn("Sequence input is empty or invalid", response.json()["detail"])

if __name__ == "__main__":
    unittest.main(verbosity=2)
