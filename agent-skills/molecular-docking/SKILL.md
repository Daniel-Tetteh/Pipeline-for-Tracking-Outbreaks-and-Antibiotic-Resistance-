# Skill: 3D Protein Folding and Chemical Compound Docking
## Capabilities
- Automate HTTP request payloads to the ESMFold API endpoint to fetch atomic coordinate files (.pdb).
- Orchestrate background command line execution blocks for molecular docking engines (like AutoDock Vina or DiffDock wrappers).
- Screen protein targets against an index list of generic antibiotic compound chemical SMILES inputs.

## Output Schema Target
Must append all docking execution calculations directly into a root dataset file named `docking_affinity_matrix.csv`.
