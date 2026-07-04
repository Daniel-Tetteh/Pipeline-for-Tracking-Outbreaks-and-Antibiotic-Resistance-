import logging
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
import sequence_engine
import structural_docking

# Setup logging for internal diagnostics
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("AMR_Pipeline")

app = FastAPI(
    title="Prokaryotic AMR Discovery Pipeline API",
    description="Backend engine for six-frame translation, ML filtering, and structural docking of AMR genes.",
    version="1.0.0"
)

class GenomeRequest(BaseModel):
    fasta_data: str

HTML_UI = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>AMR Discovery Pipeline</title>
  <meta name="description" content="Prokaryotic Antimicrobial Resistance gene discovery pipeline — six-frame translation, ML filtering, and structural docking."/>
  <link rel="preconnect" href="https://fonts.googleapis.com"/>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet"/>
  <style>
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

    :root {
      --bg: #080c14;
      --surface: #0e1420;
      --surface2: #151d2e;
      --border: #1e2d45;
      --accent: #00d4ff;
      --accent2: #7c3aed;
      --accent3: #10b981;
      --danger: #ef4444;
      --warn: #f59e0b;
      --text: #e2e8f0;
      --muted: #64748b;
      --glow: rgba(0, 212, 255, 0.15);
    }

    html { scroll-behavior: smooth; }

    body {
      background: var(--bg);
      color: var(--text);
      font-family: 'Inter', sans-serif;
      min-height: 100vh;
      overflow-x: hidden;
    }

    /* Animated background grid */
    body::before {
      content: '';
      position: fixed;
      inset: 0;
      background-image:
        linear-gradient(rgba(0,212,255,0.03) 1px, transparent 1px),
        linear-gradient(90deg, rgba(0,212,255,0.03) 1px, transparent 1px);
      background-size: 40px 40px;
      pointer-events: none;
      z-index: 0;
    }

    /* Radial glow blobs */
    .bg-blob {
      position: fixed;
      border-radius: 50%;
      filter: blur(120px);
      pointer-events: none;
      z-index: 0;
    }
    .blob1 { width: 600px; height: 600px; background: rgba(0,212,255,0.05); top: -200px; right: -100px; }
    .blob2 { width: 500px; height: 500px; background: rgba(124,58,237,0.06); bottom: -150px; left: -100px; }

    .container {
      position: relative;
      z-index: 1;
      max-width: 1100px;
      margin: 0 auto;
      padding: 0 24px 80px;
    }

    /* Header */
    header {
      display: flex;
      align-items: center;
      gap: 16px;
      padding: 40px 0 48px;
      border-bottom: 1px solid var(--border);
      margin-bottom: 48px;
    }

    .logo-icon {
      width: 48px; height: 48px;
      background: linear-gradient(135deg, var(--accent), var(--accent2));
      border-radius: 14px;
      display: flex; align-items: center; justify-content: center;
      font-size: 24px;
      box-shadow: 0 0 30px rgba(0,212,255,0.3);
      flex-shrink: 0;
    }

    .header-text h1 {
      font-size: 1.6rem;
      font-weight: 700;
      background: linear-gradient(135deg, #fff 30%, var(--accent));
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      line-height: 1.2;
    }
    .header-text p {
      color: var(--muted);
      font-size: 0.875rem;
      margin-top: 4px;
    }

    .badge {
      margin-left: auto;
      background: rgba(0,212,255,0.1);
      border: 1px solid rgba(0,212,255,0.25);
      color: var(--accent);
      font-size: 0.7rem;
      font-weight: 600;
      padding: 4px 10px;
      border-radius: 20px;
      letter-spacing: 0.05em;
      text-transform: uppercase;
    }

    /* Pipeline steps */
    .steps {
      display: grid;
      grid-template-columns: repeat(4, 1fr);
      gap: 12px;
      margin-bottom: 40px;
    }

    .step {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 16px;
      transition: border-color 0.3s, box-shadow 0.3s;
    }
    .step:hover { border-color: rgba(0,212,255,0.3); box-shadow: 0 0 20px rgba(0,212,255,0.05); }
    .step-num {
      font-size: 0.65rem;
      font-weight: 700;
      letter-spacing: 0.1em;
      color: var(--accent);
      text-transform: uppercase;
      margin-bottom: 6px;
    }
    .step-title { font-size: 0.8rem; font-weight: 600; color: var(--text); }
    .step-desc { font-size: 0.72rem; color: var(--muted); margin-top: 4px; line-height: 1.4; }

    /* Input card */
    .input-card {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 16px;
      padding: 32px;
      margin-bottom: 32px;
      transition: border-color 0.3s;
    }
    .input-card:focus-within {
      border-color: rgba(0,212,255,0.4);
      box-shadow: 0 0 40px rgba(0,212,255,0.06);
    }

    .input-label {
      display: flex;
      align-items: center;
      gap: 8px;
      font-size: 0.85rem;
      font-weight: 600;
      color: var(--text);
      margin-bottom: 12px;
    }
    .input-label span { color: var(--muted); font-weight: 400; font-size: 0.78rem; }

    #sequence-input {
      width: 100%;
      min-height: 180px;
      background: var(--surface2);
      border: 1px solid var(--border);
      border-radius: 10px;
      color: var(--text);
      font-family: 'JetBrains Mono', monospace;
      font-size: 0.8rem;
      line-height: 1.6;
      padding: 16px;
      resize: vertical;
      outline: none;
      transition: border-color 0.2s;
    }
    #sequence-input:focus { border-color: rgba(0,212,255,0.4); }
    #sequence-input::placeholder { color: var(--muted); }

    .input-footer {
      display: flex;
      align-items: center;
      justify-content: space-between;
      margin-top: 16px;
      flex-wrap: wrap;
      gap: 12px;
    }

    .char-count { font-size: 0.75rem; color: var(--muted); font-family: 'JetBrains Mono', monospace; }

    .btn-group { display: flex; gap: 10px; }

    .btn {
      display: flex; align-items: center; gap: 8px;
      padding: 10px 22px;
      border-radius: 10px;
      font-size: 0.85rem;
      font-weight: 600;
      cursor: pointer;
      border: none;
      transition: all 0.2s;
      font-family: 'Inter', sans-serif;
    }
    .btn-primary {
      background: linear-gradient(135deg, var(--accent), #0099bb);
      color: #000;
      box-shadow: 0 0 20px rgba(0,212,255,0.3);
    }
    .btn-primary:hover { transform: translateY(-1px); box-shadow: 0 4px 30px rgba(0,212,255,0.45); }
    .btn-primary:active { transform: translateY(0); }
    .btn-primary:disabled { opacity: 0.5; cursor: not-allowed; transform: none; }

    .btn-ghost {
      background: transparent;
      color: var(--muted);
      border: 1px solid var(--border);
    }
    .btn-ghost:hover { color: var(--text); border-color: rgba(255,255,255,0.2); }

    /* Status bar */
    #status-bar {
      display: none;
      align-items: center;
      gap: 12px;
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 14px 20px;
      margin-bottom: 24px;
      font-size: 0.82rem;
    }

    .spinner {
      width: 18px; height: 18px;
      border: 2px solid rgba(0,212,255,0.2);
      border-top-color: var(--accent);
      border-radius: 50%;
      animation: spin 0.8s linear infinite;
      flex-shrink: 0;
    }
    @keyframes spin { to { transform: rotate(360deg); } }

    .pulse-dot {
      width: 10px; height: 10px;
      border-radius: 50%;
      background: var(--accent3);
      animation: pulse 1.2s ease-in-out infinite;
      flex-shrink: 0;
    }
    @keyframes pulse { 0%, 100% { opacity: 1; transform: scale(1); } 50% { opacity: 0.5; transform: scale(0.7); } }

    /* Results */
    #results-container { display: none; }

    .results-header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      margin-bottom: 20px;
    }
    .results-title {
      font-size: 1rem;
      font-weight: 600;
    }
    .results-meta { font-size: 0.78rem; color: var(--muted); }

    .orf-card {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 14px;
      margin-bottom: 20px;
      overflow: hidden;
      transition: border-color 0.2s;
      animation: fadeIn 0.4s ease forwards;
      opacity: 0;
    }
    .orf-card:hover { border-color: rgba(0,212,255,0.25); }
    @keyframes fadeIn { to { opacity: 1; } }

    .orf-header {
      display: flex;
      align-items: center;
      gap: 12px;
      padding: 18px 24px;
      border-bottom: 1px solid var(--border);
      cursor: pointer;
      user-select: none;
    }

    .status-pill {
      font-size: 0.68rem;
      font-weight: 700;
      letter-spacing: 0.06em;
      text-transform: uppercase;
      padding: 3px 10px;
      border-radius: 20px;
    }
    .pill-mutated { background: rgba(245,158,11,0.15); color: var(--warn); border: 1px solid rgba(245,158,11,0.3); }
    .pill-novel   { background: rgba(16,185,129,0.15); color: var(--accent3); border: 1px solid rgba(16,185,129,0.3); }
    .pill-pseudo  { background: rgba(239,68,68,0.12); color: var(--danger); border: 1px solid rgba(239,68,68,0.25); }

    .orf-title { font-size: 0.88rem; font-weight: 600; flex: 1; }
    .orf-meta  { font-size: 0.75rem; color: var(--muted); }
    .chevron { color: var(--muted); transition: transform 0.2s; font-size: 0.9rem; }
    .chevron.open { transform: rotate(180deg); }

    .orf-body { padding: 20px 24px; display: none; }
    .orf-body.visible { display: block; }

    .info-grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
      gap: 12px;
      margin-bottom: 20px;
    }
    .info-cell {
      background: var(--surface2);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px 14px;
    }
    .info-cell-label { font-size: 0.68rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }
    .info-cell-value { font-size: 0.85rem; font-weight: 600; font-family: 'JetBrains Mono', monospace; }

    .aa-seq-box {
      background: var(--surface2);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 14px;
      font-family: 'JetBrains Mono', monospace;
      font-size: 0.75rem;
      line-height: 1.8;
      word-break: break-all;
      color: #94d2ff;
      margin-bottom: 20px;
      max-height: 100px;
      overflow-y: auto;
    }

    .docking-title {
      font-size: 0.8rem;
      font-weight: 600;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.05em;
      margin-bottom: 12px;
    }

    .docking-grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(160px, 1fr));
      gap: 10px;
    }

    .docking-cell {
      background: var(--surface2);
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 12px 14px;
      transition: border-color 0.2s, transform 0.2s;
    }
    .docking-cell:hover { transform: translateY(-2px); border-color: rgba(0,212,255,0.25); }
    .docking-compound { font-size: 0.72rem; color: var(--muted); margin-bottom: 4px; }
    .docking-affinity {
      font-size: 1rem;
      font-weight: 700;
      font-family: 'JetBrains Mono', monospace;
    }
    .affinity-high   { color: var(--danger); }
    .affinity-medium { color: var(--warn); }
    .affinity-low    { color: var(--accent3); }
    .docking-class { font-size: 0.65rem; color: var(--muted); margin-top: 2px; }

    /* Error card */
    .error-card {
      background: rgba(239,68,68,0.08);
      border: 1px solid rgba(239,68,68,0.3);
      border-radius: 12px;
      padding: 20px 24px;
      color: #fca5a5;
      font-size: 0.85rem;
      display: flex; gap: 12px; align-items: flex-start;
      animation: fadeIn 0.3s ease forwards;
    }
    .error-icon { font-size: 1.2rem; flex-shrink: 0; }

    /* Empty state */
    .empty-state {
      text-align: center;
      padding: 60px 20px;
      color: var(--muted);
    }
    .empty-state-icon { font-size: 3rem; margin-bottom: 16px; opacity: 0.5; }
    .empty-state h3 { font-size: 1rem; font-weight: 600; color: var(--text); margin-bottom: 8px; }
    .empty-state p { font-size: 0.82rem; }

    /* Scrollbar */
    ::-webkit-scrollbar { width: 6px; height: 6px; }
    ::-webkit-scrollbar-track { background: var(--surface2); }
    ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }

    @media (max-width: 768px) {
      .steps { grid-template-columns: repeat(2, 1fr); }
      .badge { display: none; }
    }
    @media (max-width: 480px) {
      .steps { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <div class="bg-blob blob1"></div>
  <div class="bg-blob blob2"></div>

  <div class="container">
    <header>
      <div class="logo-icon">🧬</div>
      <div class="header-text">
        <h1>AMR Discovery Pipeline</h1>
        <p>Prokaryotic Antimicrobial Resistance Gene Analysis</p>
      </div>
      <div class="badge">v1.0.0 Live</div>
    </header>

    <!-- Pipeline steps overview -->
    <div class="steps">
      <div class="step">
        <div class="step-num">Step 01</div>
        <div class="step-title">IUPAC Sanitization</div>
        <div class="step-desc">Validates characters against A, T, G, C, N codes and aborts on injection.</div>
      </div>
      <div class="step">
        <div class="step-num">Step 02</div>
        <div class="step-title">6-Frame Translation</div>
        <div class="step-desc">Translates all forward and reverse frames, extracts ORFs ≥ 30 aa.</div>
      </div>
      <div class="step">
        <div class="step-num">Step 03</div>
        <div class="step-title">ML Filter + Homology</div>
        <div class="step-desc">Random Forest noise removal and Needleman-Wunsch AMR alignment.</div>
      </div>
      <div class="step">
        <div class="step-num">Step 04</div>
        <div class="step-title">Structural Docking</div>
        <div class="step-desc">ESMFold structure prediction and antibiotic compound screening.</div>
      </div>
    </div>

    <!-- Input card -->
    <div class="input-card">
      <div class="input-label">
        Paste DNA Sequence
        <span>— raw nucleotides or FASTA format (only A, T, G, C, N allowed)</span>
      </div>
      <textarea
        id="sequence-input"
        spellcheck="false"
        autocomplete="off"
        placeholder=">Sample_Isolate&#10;ATGTCTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCT&#10;CACCCAGAAACGCTGGTAAAGTAAAAGATGCGGAAGATCAGTTGGGTGCACGAGTGGG..."
      ></textarea>
      <div class="input-footer">
        <span class="char-count" id="char-count">0 nucleotides</span>
        <div class="btn-group">
          <button class="btn btn-ghost" id="clear-btn" onclick="clearInput()">✕ Clear</button>
          <button class="btn btn-primary" id="analyze-btn" onclick="analyzeSequence()">
            <span id="btn-icon">⚡</span>
            <span id="btn-text">Run Analysis</span>
          </button>
        </div>
      </div>
    </div>

    <!-- Status bar -->
    <div id="status-bar">
      <div class="spinner" id="status-spinner"></div>
      <div class="pulse-dot" id="status-dot" style="display:none"></div>
      <span id="status-text">Initializing pipeline...</span>
    </div>

    <!-- Results -->
    <div id="results-container">
      <div class="results-header">
        <div class="results-title">Analysis Results</div>
        <div class="results-meta" id="results-meta"></div>
      </div>
      <div id="results-body"></div>
    </div>
  </div>

  <script>
    const input = document.getElementById('sequence-input');
    const charCount = document.getElementById('char-count');
    const analyzeBtn = document.getElementById('analyze-btn');
    const btnText = document.getElementById('btn-text');
    const btnIcon = document.getElementById('btn-icon');
    const statusBar = document.getElementById('status-bar');
    const statusText = document.getElementById('status-text');
    const statusSpinner = document.getElementById('status-spinner');
    const statusDot = document.getElementById('status-dot');
    const resultsContainer = document.getElementById('results-container');
    const resultsMeta = document.getElementById('results-meta');
    const resultsBody = document.getElementById('results-body');

    // Live nucleotide counter
    input.addEventListener('input', () => {
      const raw = input.value.replace(/>.*?\\n/g, '').replace(/\\s/g, '');
      charCount.textContent = raw.length.toLocaleString() + ' nucleotides';
    });

    function clearInput() {
      input.value = '';
      charCount.textContent = '0 nucleotides';
      resultsContainer.style.display = 'none';
      statusBar.style.display = 'none';
      resultsBody.innerHTML = '';
    }

    function setStatus(msg, loading = true) {
      statusBar.style.display = 'flex';
      statusText.textContent = msg;
      statusSpinner.style.display = loading ? 'block' : 'none';
      statusDot.style.display = loading ? 'none' : 'block';
    }

    function affinityClass(val) {
      if (val <= -7.5) return 'affinity-high';
      if (val <= -5.5) return 'affinity-medium';
      return 'affinity-low';
    }

    function pillClass(orf) {
      if (orf.is_pseudogene) return 'pill-pseudo';
      if (orf.status === 'MUTATED_VARIANT') return 'pill-mutated';
      return 'pill-novel';
    }
    function pillLabel(orf) {
      if (orf.is_pseudogene) return '☠ Pseudogene';
      if (orf.status === 'MUTATED_VARIANT') return '⚠ Mutated Variant';
      return '✦ Novel Protein';
    }

    function compoundClass(name) {
      const map = {
        'Penicillin G': 'beta-lactam', 'Ampicillin': 'beta-lactam',
        'Ciprofloxacin': 'fluoroquinolone', 'Tetracycline': 'tetracycline',
        'Erythromycin': 'macrolide', 'Gentamicin': 'aminoglycoside'
      };
      return map[name] || '';
    }

    function toggleOrfBody(index) {
      const body = document.getElementById('orf-body-' + index);
      const chevron = document.getElementById('chevron-' + index);
      body.classList.toggle('visible');
      chevron.classList.toggle('open');
    }

    function renderResults(data) {
      resultsContainer.style.display = 'block';
      const count = data.results.length;
      resultsMeta.textContent = count + ' coding sequence' + (count !== 1 ? 's' : '') + ' identified';

      if (count === 0) {
        resultsBody.innerHTML = `
          <div class="empty-state">
            <div class="empty-state-icon">🔬</div>
            <h3>No ORFs Detected</h3>
            <p>No coding sequences above 30 amino acids were found in this input. Try a longer sequence.</p>
          </div>`;
        return;
      }

      resultsBody.innerHTML = data.results.map((orf, i) => {
        const dockingHTML = orf.docking_results.map(d => `
          <div class="docking-cell">
            <div class="docking-compound">${d.compound_name}</div>
            <div class="docking-affinity ${affinityClass(d.binding_affinity_kcal_mol)}">${d.binding_affinity_kcal_mol} kcal/mol</div>
            <div class="docking-class">${compoundClass(d.compound_name)}</div>
          </div>`).join('');

        const frameLabel = (orf.is_reverse ? 'Rev ' : 'Fwd ') + 'Frame ' + (orf.frame + 1);
        const identPct = (orf.sequence_identity * 100).toFixed(1) + '%';
        const refLabel = orf.reference_match || '—';

        const delay = i * 80;

        return `
          <div class="orf-card" style="animation-delay:${delay}ms">
            <div class="orf-header" onclick="toggleOrfBody(${i})">
              <div class="status-pill ${pillClass(orf)}">${pillLabel(orf)}</div>
              <div class="orf-title">ORF ${i + 1} &mdash; ${orf.amino_acid_sequence.length} aa &mdash; ${frameLabel}</div>
              <div class="orf-meta">${refLabel}</div>
              <div class="chevron" id="chevron-${i}">▼</div>
            </div>
            <div class="orf-body" id="orf-body-${i}">
              <div class="info-grid">
                <div class="info-cell">
                  <div class="info-cell-label">Status</div>
                  <div class="info-cell-value">${orf.status}</div>
                </div>
                <div class="info-cell">
                  <div class="info-cell-label">Reference Match</div>
                  <div class="info-cell-value">${refLabel}</div>
                </div>
                <div class="info-cell">
                  <div class="info-cell-label">Sequence Identity</div>
                  <div class="info-cell-value">${identPct}</div>
                </div>
                <div class="info-cell">
                  <div class="info-cell-label">Reading Frame</div>
                  <div class="info-cell-value">${frameLabel}</div>
                </div>
                <div class="info-cell">
                  <div class="info-cell-label">Action Required</div>
                  <div class="info-cell-value" style="color:${orf.action_required ? '#f59e0b' : '#10b981'}">${orf.action_required ? 'YES' : 'NO'}</div>
                </div>
                <div class="info-cell">
                  <div class="info-cell-label">Pseudogene</div>
                  <div class="info-cell-value" style="color:${orf.is_pseudogene ? '#ef4444' : '#10b981'}">${orf.is_pseudogene ? 'YES' : 'NO'}</div>
                </div>
              </div>
              <div class="aa-seq-box">${orf.amino_acid_sequence}</div>
              <div class="docking-title">Antibiotic Binding Affinities</div>
              <div class="docking-grid">${dockingHTML}</div>
              ${orf.notes ? '<div style="margin-top:14px;font-size:0.78rem;color:#f59e0b;background:rgba(245,158,11,0.08);border:1px solid rgba(245,158,11,0.2);border-radius:8px;padding:10px 14px;">⚠ ' + orf.notes + '</div>' : ''}
            </div>
          </div>`;
      }).join('');

      // Auto-expand first orf
      if (count > 0) toggleOrfBody(0);
    }

    async function analyzeSequence() {
      const seq = input.value.trim();
      if (!seq) return;

      analyzeBtn.disabled = true;
      btnIcon.textContent = '⏳';
      btnText.textContent = 'Analyzing...';
      resultsContainer.style.display = 'none';
      resultsBody.innerHTML = '';

      const steps = [
        'Sanitizing IUPAC nucleotide input...',
        'Running 6-frame translation...',
        'Applying ML noise filter...',
        'Performing homology alignment...',
        'Predicting protein structure (ESMFold)...',
        'Screening antibiotic compound docking...',
        'Compiling results...'
      ];
      let stepIdx = 0;
      setStatus(steps[0]);
      const stepTimer = setInterval(() => {
        stepIdx = Math.min(stepIdx + 1, steps.length - 1);
        setStatus(steps[stepIdx]);
      }, 800);

      try {
        const res = await fetch('/api/v1/analyze_genome', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ fasta_data: seq })
        });
        clearInterval(stepTimer);
        const data = await res.json();

        if (!res.ok) {
          statusBar.style.display = 'none';
          resultsContainer.style.display = 'block';
          resultsMeta.textContent = '';
          resultsBody.innerHTML = `
            <div class="error-card">
              <div class="error-icon">🛡</div>
              <div><strong>Pipeline Rejected Input</strong><br/>${data.detail}</div>
            </div>`;
        } else {
          setStatus('Analysis complete — ' + data.results.length + ' proteins processed.', false);
          renderResults(data);
        }
      } catch (err) {
        clearInterval(stepTimer);
        statusBar.style.display = 'none';
        resultsContainer.style.display = 'block';
        resultsMeta.textContent = '';
        resultsBody.innerHTML = `
          <div class="error-card">
            <div class="error-icon">⚡</div>
            <div><strong>Network Error</strong><br/>Could not reach the pipeline server. Is uvicorn running?</div>
          </div>`;
      } finally {
        analyzeBtn.disabled = false;
        btnIcon.textContent = '⚡';
        btnText.textContent = 'Run Analysis';
      }
    }

    // Allow Ctrl+Enter to submit
    input.addEventListener('keydown', e => {
      if (e.ctrlKey && e.key === 'Enter') analyzeSequence();
    });
  </script>
</body>
</html>"""


@app.get("/", response_class=HTMLResponse)
async def serve_ui():
    """Serves the interactive AMR Pipeline web UI."""
    return HTML_UI


@app.post("/api/v1/analyze_genome")
async def analyze_genome(request: GenomeRequest):
    """
    Exposes a single POST route accepting raw genome fasta configurations,
    running them through IUPAC validation, translation, ML filtration,
    homology checks, and molecular docking.
    """
    try:
        raw_sequence = request.fasta_data
        logger.info(f"Received genome sequence of length {len(raw_sequence)}")

        analyzed_orfs = sequence_engine.extract_and_analyze_orfs(raw_sequence)
        logger.info(f"Extracted and filtered {len(analyzed_orfs)} biological coding sequences")

        results = []

        for orf in analyzed_orfs:
            aa_seq = orf["amino_acid_sequence"]
            ref_match = orf["reference_match"]
            seq_identity = orf["sequence_identity"]

            logger.info(f"Screening target sequence ({len(aa_seq)} aa) against antibiotic compounds index...")

            docking_scores = structural_docking.screen_compounds(
                amino_acid_sequence=aa_seq,
                ref_match=ref_match,
                seq_identity=seq_identity
            )

            orf_response = {
                "amino_acid_sequence": aa_seq,
                "status": orf["status"],
                "action_required": orf["action_required"],
                "sequence_identity": seq_identity,
                "reference_match": ref_match,
                "frame": orf["frame"],
                "is_reverse": orf["is_reverse"],
                "is_pseudogene": orf["is_pseudogene"],
                "docking_results": docking_scores
            }
            if "notes" in orf:
                orf_response["notes"] = orf["notes"]

            results.append(orf_response)

        return {
            "success": True,
            "message": f"Processed {len(results)} target proteins.",
            "results": results
        }

    except ValueError as e:
        logger.warning(f"Sanitization alert: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))

    except Exception as e:
        logger.error(f"Internal Pipeline Error: {str(e)}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail="An internal computational pipeline error occurred. The system trace was sanitized for safety."
        )
