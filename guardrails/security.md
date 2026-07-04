# Agent Security & Quality Guardrails
## STRIDE Threat Mitigation Rules
- **Tampering / Injection**: Sanitize incoming sequence variables. If strings contain elements outside valid IUPAC nucleotide codes (A, T, G, C, N), execution must abort.
- **Information Disclosure**: Catch script runtime tracebacks using clean try/except blocks to stop internal system layout leaks.
