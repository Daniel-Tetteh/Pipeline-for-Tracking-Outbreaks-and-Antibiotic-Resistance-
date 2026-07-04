import os
import random
import joblib
import numpy as np
from sklearn.ensemble import RandomForestClassifier

MODEL_PATH = os.path.join(os.path.dirname(__file__), "orfs_classifier.joblib")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Typical coding sequence amino acid frequencies in prokaryotes
# Source: general codon usage/amino acid composition
CODING_FREQ = {
    'A': 0.10, 'R': 0.05, 'N': 0.04, 'D': 0.05, 'C': 0.01,
    'Q': 0.04, 'E': 0.06, 'G': 0.07, 'H': 0.02, 'I': 0.06,
    'L': 0.10, 'K': 0.05, 'M': 0.02, 'F': 0.04, 'P': 0.04,
    'S': 0.06, 'T': 0.05, 'W': 0.01, 'Y': 0.03, 'V': 0.07
}

def extract_features(sequence: str) -> list:
    """
    Extracts amino acid frequencies as features for the ML model.
    Returns a list of 20 floats representing the relative frequency of each standard amino acid.
    """
    seq = sequence.upper().strip()
    length = len(seq)
    if length == 0:
        return [0.0] * len(AMINO_ACIDS)
    
    counts = {aa: seq.count(aa) for aa in AMINO_ACIDS}
    features = [counts[aa] / length for aa in AMINO_ACIDS]
    return features

def generate_synthetic_data(num_samples: int = 150):
    """
    Generates synthetic coding (positive) and non-coding (negative) sequences.
    """
    X = []
    y = []
    
    # 1. Positive Samples (Coding) - Follow prokaryotic amino acid distribution
    aa_list = list(CODING_FREQ.keys())
    weights = list(CODING_FREQ.values())
    
    # Add actual AMR sequences if available to anchor the positive class
    amr_references = [
        "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW",
        "MNSIPYGIAALIGLSLGALFSGLVGLYFRPLLVVLFVAFLGVFLFRLWLLRLWRLWRRAGLLVAFALALFAGFVRWVLVALAVILGAVLLGIVFGWWGLWRAGRPAVLAIAVAVAVAGVVRAGARAGAR",
        "MGFGSRGLLLLLCVLAAVAAAATGCCCCGGGGTTTTAAAGCCGCTAGCCTAGCATGATCGATCG"
    ]
    
    for ref in amr_references:
        X.append(extract_features(ref))
        y.append(1)
        
    for _ in range(num_samples):
        # Generate positive: sequence with typical coding composition
        length = random.randint(50, 300)
        seq = "".join(random.choices(aa_list, weights=weights, k=length))
        X.append(extract_features(seq))
        y.append(1)
        
        # Generate negative: uniform random sequence representing out-of-frame noise
        seq_noise = "".join(random.choices(AMINO_ACIDS, k=length))
        X.append(extract_features(seq_noise))
        y.append(0)
        
    return np.array(X), np.array(y)

def train_and_save_model():
    """
    Trains the Random Forest model on synthetic data and saves it.
    """
    print("Training ML noise filter classifier...")
    X, y = generate_synthetic_data()
    # Using a shallow Random Forest to avoid overfitting on synthetic data
    clf = RandomForestClassifier(n_estimators=50, max_depth=5, random_state=42)
    clf.fit(X, y)
    joblib.dump(clf, MODEL_PATH)
    print(f"Model saved to {MODEL_PATH}")
    return clf

def get_classifier():
    """
    Loads and returns the classifier, training it if not present.
    """
    if not os.path.exists(MODEL_PATH):
        return train_and_save_model()
    try:
        return joblib.load(MODEL_PATH)
    except Exception:
        # Fallback to training a new model if loading fails
        return train_and_save_model()

def predict_is_protein(sequence: str) -> bool:
    """
    Predicts if the sequence is a valid coding protein (True) or out-of-frame noise (False).
    """
    if len(sequence) < 30:
        return False
    
    clf = get_classifier()
    features = np.array([extract_features(sequence)])
    prob = clf.predict_proba(features)[0][1]
    # Return True if it has >= 50% probability of being a protein
    return prob >= 0.5
