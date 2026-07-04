# Skill: Six-Frame Translation and Genomic Alignment Filter
## Capabilities
- Read raw prokaryotic DNA sequence inputs.
- Translate all 6 forward and reverse reading frames utilizing BioPython.
- Filter biological noise by discarding sequences under 30 amino acids.
- Perform high-speed sequence comparison checks against local reference databases.

## Output Schema Target
Must return a clean JSON sequence block containing:
- `amino_acid_sequence`: string
- `status`: "MUTATED_VARIANT" or "NOVEL_PROTEIN"
- `action_required`: true/false
