# AMR-Agent: Autonomous Six-Frame Variant Proteome Pipeline for Tracking Outbreaks and Antibiotic Resistance


## 🧬 Project Description

Antimicrobial Resistance (AMR) and the rapid emergence of novel bacterial pathogens represent one of the most critical challenges to global public health. Traditional prokaryotic gene-prediction tools rely on rigid, heuristic algorithms that search for standard, clean Open Reading Frames (ORFs). However, when pathogens mutate to survive antibiotic exposure, their genomes often undergo frameshifts, alternative start-codon selections, and pseudogenization. Standard workflows routinely drop these highly critical variant sequences as "non-coding junk noise," leaving a dangerous blind spot in clinical biodefense.

**AMR-Agent** is an entirely autonomous, spec-driven computational pipeline built inside the **Google Antigravity 2.0 IDE** framework [comp]. Operating on the principles of file-based context engineering and the Model Context Protocol (MCP), the agent acts as an automated genomic investigator to bridge the gap between raw genomic sequence data and structural drug interaction analytics.

### How the Pipeline Works:
1. **6-Frame Variant Isolation:** The pipeline reads raw, unannotated prokaryotic DNA fasta records and performs a full six-frame window translation (3 forward frames, 3 reverse frames) to capture all theoretical translation possibilities.
2. **Reference Proteome Subtraction:** It compares the output matrix against a high-fidelity reference database (such as *Escherichia coli* K-12 MG1655). Identical sequences are logged and terminated, while novel inserts and mutated partial matches are pushed forward.
3. **Machine Learning Noise Filtering:** It extracts structural sequence metrics (amino acid distribution, length, GC content) and processes them through an integrated XGBoost classifier to distinguish functional sequences from random frame noise. It embeds specific heuristics to flag and preserve truncated pseudogenes caused by frameshifts.
4. **Structural AI & High-Throughput Docking:** Approved novel or variant target proteins are automatically dispatched to the ESMFold API to generate 3D atomic structures (`.pdb` coordinate files). The agent then runs an automated command-line subprocess to dock these structures against a chemical library of generic antibiotics (using SMILES strings) to output binding affinity matrices (`docking_affinity_matrix.csv`).

By shifting security left and embedding automated STRIDE threat mitigation boundaries, this project demonstrates how an autonomous software workforce can compress pathogen evaluation timelines from weeks of wet-lab cultivation down to under a single minute of agentic computation.

---

## 🏗️ Repository Architecture

This workspace utilizes the standard file-based skill and constraint configuration blueprint expected by the Google Antigravity core engine [comp]:

```text
amr-antigravity-pipeline/
├── agent-skills/
│   ├── extract-orfs/
│   │   └── SKILL.md          <- Day 3 Context Engineering: Translation & Alignment
│   └── molecular-docking/
│       └── SKILL.md          <- Day 3 Context Engineering: Folding & Docking Loops
├── guardrails/
│   └── security.md           <- Day 4 STRIDE Threat Mitigation Boundaries
├── specs/
│   └── production_spec.md    <- Day 5 Production App Service Blueprint
├── data/
│   └── test_genome.fasta     <- Sample Isolate Mock Input DNA Sequence
├── sequence_engine.py        <- Core Six-Frame Translation Module (Agent Generated)
├── ml_filter.py              <- Noise Filtering & Pseudogene Logic (Agent Generated)
├── structural_docking.py     <- ESMFold API & Docking Aggregation (Agent Generated)
├── app.py                    <- Production FastAPI Gateway Wrapper (Agent Generated)
├── requirements.txt          <- Python Sandbox Dependency Manifest
└── README.md                 <- System Setup & Documentation Hub
```

---
## 🛠️ Local Installation & Setup Instructions

Execute these terminal instructions inside your workspace to reproduce the sandbox environment and run the pipeline locally:

### 1. Clone the Repository
```bash
git clone https://github.comYOUR_GITHUB_USERNAME/amr-antigravity-pipeline.git
cd amr-antigravity-pipeline
```

### 2. Configure Your Virtual Environment
Create a virtual python isolation sandbox and install the system dependencies:
```bash
python -m venv amr_env
source amr_env/bin/activate  # On Windows: amr_env\Scripts\activate
pip install -r requirements.txt
```

### 3. Run the Sandbox Simulation Pipeline
Run the entry point application to execute the full translation, filtering, and docking analysis loop against the provided sample fasta genome file:
```bash
python app.py
```
