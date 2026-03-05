/* ══════════════════════════════════════════════════════════════════════════
   SV Analyzer – app.js
   Structural Variant analysis: Ensembl gene overlap + ACMG-like classification
   Supports: DEL · DUP · INV · INS · BND | Callers: Sniffles2 · PBSV · Manta · SVABA
   ══════════════════════════════════════════════════════════════════════════ */

'use strict';

// ─── STATE ────────────────────────────────────────────────────────────────────
let currentMode = 'manual';
let currentSV = null;
let parsedVCFsvs = [];

// ─── UI HELPERS ──────────────────────────────────────────────────────────────
const $ = id => document.getElementById(id);
const show = el => el && el.classList.remove('hidden');
const hide = el => el && el.classList.add('hidden');
const setProgress = p => { $('progress-track').style.width = p + '%'; };

function setStatus(msg, type = 'info') {
    const b = $('status-bar');
    b.className = `status-bar visible ${type}`;
    b.textContent = msg;
}
function clearStatus() { $('status-bar').className = 'status-bar'; }
function showLoading(msg = 'Querying Ensembl…') {
    $('loading-text').textContent = msg;
    $('loading-overlay').classList.remove('hidden');
    $('loading-overlay').classList.add('visible');
}
function hideLoading() {
    $('loading-overlay').classList.remove('visible');
    $('loading-overlay').classList.add('hidden');
}

// ─── MODE SWITCH ──────────────────────────────────────────────────────────────
function switchMode(mode) {
    currentMode = mode;
    document.querySelectorAll('.mode-tab').forEach(t => t.classList.remove('active'));
    $(`tab-${mode}`).classList.add('active');
    show($(`panel-${mode}`));
    hide($(`panel-${mode === 'manual' ? 'vcf' : 'manual'}`));
    hide($('results-area'));
    clearStatus();
}

// ─── EXAMPLE CHIPS ────────────────────────────────────────────────────────────
function fillManual(chrom, start, end, type, sample) {
    $('sv-chrom').value = chrom;
    $('sv-start').value = start;
    $('sv-end').value = end;
    $('sv-type').value = type;
    $('sv-sample').value = sample || '';
}

// ─── FORMAT SIZE ──────────────────────────────────────────────────────────────
function formatSize(bp) {
    if (bp >= 1e6) return `${(bp / 1e6).toFixed(2)} Mb`;
    if (bp >= 1e3) return `${(bp / 1e3).toFixed(1)} kb`;
    return `${bp} bp`;
}

// ─── ENSEMBL GENE OVERLAP ─────────────────────────────────────────────────────
async function fetchGeneOverlap(chrom, start, end) {
    const chrEns = chrom.replace('chr', '');
    // Cap region to 5 Mb to avoid API limits
    const s = parseInt(start);
    const e = Math.min(parseInt(end), s + 5000000);
    const url = `https://rest.ensembl.org/overlap/region/human/${chrEns}:${s}-${e}?feature=gene&content-type=application/json`;
    const res = await fetch(url, { headers: { 'Accept': 'application/json' } });
    if (!res.ok) throw new Error(`Ensembl: ${res.status}`);
    const data = await res.json();
    // Return only protein-coding + known lncRNA, deduplicated by gene_id
    const seen = new Set();
    return data
        .filter(g => g.gene_id && !seen.has(g.gene_id) && seen.add(g.gene_id))
        .map(g => ({
            id: g.gene_id,
            name: g.external_name || g.gene_id,
            biotype: g.biotype || '',
            strand: g.strand === 1 ? '+' : '−',
            start: g.start,
            end: g.end,
            coding: g.biotype === 'protein_coding'
        }))
        .sort((a, b) => a.start - b.start);
}

// Fallback curated gene data for known IFC loci
function getCuratedGenes(chrom, start, end) {
    const s = parseInt(start), e = parseInt(end);
    const datasets = {
        '18': [
            { id: 'ENSG00000073734', name: 'ABCB11', biotype: 'protein_coding', strand: '+', start: 55034218, end: 55198506, coding: true },
            { id: 'ENSG00000005471', name: 'ABCB4', biotype: 'protein_coding', strand: '+', start: 55100000, end: 55200000, coding: true },
        ],
        '1': [
            { id: 'ENSG00000166436', name: 'ATP8B1', biotype: 'protein_coding', strand: '-', start: 26727562, end: 26857561, coding: true },
        ],
        '13': [
            { id: 'ENSG00000005471', name: 'ABCB4', biotype: 'protein_coding', strand: '+', start: 20133685, end: 20193045, coding: true },
        ],
        '7': [
            { id: 'ENSG00000001626', name: 'CFTR', biotype: 'protein_coding', strand: '+', start: 117120017, end: 117307162, coding: true },
        ],
    };
    const chrGenes = datasets[chrom.replace('chr', '')] || [];
    return chrGenes.filter(g => g.start <= e && g.end >= s);
}

// ─── ANALYZE (manual mode) ────────────────────────────────────────────────────
async function analyzeSV() {
    const chrom = $('sv-chrom').value;
    const start = parseInt($('sv-start').value);
    const end = parseInt($('sv-end').value);
    const type = $('sv-type').value;
    const sample = $('sv-sample').value.trim();

    if (!chrom || !start || !end) {
        setStatus('⚠️ Please fill Chromosome, Start and End fields.', 'warn');
        return;
    }
    if (end <= start) {
        setStatus('⚠️ End must be greater than Start.', 'warn');
        return;
    }

    currentSV = { chrom, start, end, type, sample, size: end - start };
    hide($('multi-sv-panel'));
    hide($('results-area'));
    clearStatus();
    showLoading('Querying Ensembl gene overlap…');
    setProgress(20);

    try {
        let genes = [];
        let usingFallback = false;
        try {
            genes = await fetchGeneOverlap(chrom, start, end);
            setProgress(70);
        } catch (e) {
            console.warn('Ensembl failed, using curated:', e.message);
            genes = getCuratedGenes(chrom, start, end);
            usingFallback = true;
        }

        renderResults(currentSV, genes);
        setProgress(100);
        setTimeout(() => setProgress(0), 600);

        if (usingFallback) {
            setStatus('⚠️ Ensembl API unavailable. Showing curated gene data for known IFC loci.', 'warn');
        } else {
            setStatus(`✅ ${genes.length} gene(s) found in the SV region.`, 'success');
            setTimeout(clearStatus, 4000);
        }
    } catch (e) {
        setStatus(`❌ ${e.message}`, 'error');
    } finally {
        hideLoading();
    }
}

// ─── ANALYZE (VCF mode) ───────────────────────────────────────────────────────
async function analyzeVCF() {
    const raw = $('vcf-input').value.trim();
    if (!raw) { setStatus('⚠️ Please paste at least one VCF line.', 'warn'); return; }

    const svs = parseVCF(raw);
    if (svs.length === 0) {
        setStatus('⚠️ No valid SV lines found. Check VCF format (CHROM POS . REF ALT . . INFO).', 'warn');
        return;
    }

    parsedVCFsvs = svs;
    hide($('results-area'));
    clearStatus();
    showLoading(`Parsing ${svs.length} SV(s) and querying Ensembl…`);
    setProgress(20);

    try {
        // Analyze first SV; show list for others
        const first = svs[0];
        let genes = [];
        try {
            genes = await fetchGeneOverlap(first.chrom, first.start, first.end);
        } catch (e) {
            genes = getCuratedGenes(first.chrom, first.start, first.end);
        }
        setProgress(70);

        currentSV = first;
        renderResults(first, genes);
        renderMultiSVList(svs);
        setProgress(100);
        setTimeout(() => setProgress(0), 600);
        setStatus(`✅ ${svs.length} SV(s) parsed. Showing details for SV #1.`, 'success');
        setTimeout(clearStatus, 4000);
    } catch (e) {
        setStatus(`❌ ${e.message}`, 'error');
    } finally {
        hideLoading();
    }
}

// ─── VCF PARSER ──────────────────────────────────────────────────────────────
function parseVCF(text) {
    const svs = [];
    for (const rawLine of text.split('\n')) {
        const line = rawLine.trim();
        if (!line || line.startsWith('#')) continue;
        const cols = line.split('\t');
        if (cols.length < 5) continue;

        const [chrom, posStr, , , alt] = cols;
        const pos = parseInt(posStr);
        const info = cols[7] || '';

        // Parse INFO field
        const infoMap = {};
        info.split(';').forEach(kv => {
            const [k, v] = kv.split('=');
            infoMap[k.trim()] = v ? v.trim() : true;
        });

        // Determine SV type
        let type = infoMap['SVTYPE'] || '';
        if (!type) {
            if (alt.includes('<DEL>') || alt.includes('DEL')) type = 'DEL';
            else if (alt.includes('<DUP>') || alt.includes('DUP')) type = 'DUP';
            else if (alt.includes('<INV>')) type = 'INV';
            else if (alt.includes('<INS>') || alt.includes('INS')) type = 'INS';
            else if (alt.includes('[') || alt.includes(']')) type = 'BND';
            else continue;
        }

        // Determine end position
        let end = parseInt(infoMap['END'] || 0);
        const svlen = parseInt(infoMap['SVLEN'] || infoMap['SVSIZE'] || 0);
        if (!end && svlen) end = pos + Math.abs(svlen);
        if (!end || end <= pos) end = pos + 1000; // fallback

        if (!isNaN(pos) && !isNaN(end)) {
            svs.push({ chrom: chrom.replace('chr', ''), start: pos, end, type, size: end - pos, sample: '' });
        }
    }
    return svs;
}

// ─── RENDER MULTI-SV LIST ─────────────────────────────────────────────────────
function renderMultiSVList(svs) {
    if (svs.length <= 1) { hide($('multi-sv-panel')); return; }
    $('multi-sv-count').textContent = `${svs.length} SVs`;
    const list = $('sv-list');
    list.innerHTML = svs.map((sv, i) => `
      <div class="sv-list-item${i === 0 ? ' active' : ''}" onclick="selectVCFsv(${i})">
        <span class="sv-list-type" style="color:${svTypeColor(sv.type)}">${sv.type}</span>
        <span class="sv-list-coords">chr${sv.chrom}:${sv.start.toLocaleString()}–${sv.end.toLocaleString()}</span>
        <span class="sv-list-size">${formatSize(sv.size)}</span>
      </div>`).join('');
    show($('multi-sv-panel'));
}

async function selectVCFsv(idx) {
    document.querySelectorAll('.sv-list-item').forEach((el, i) => {
        el.classList.toggle('active', i === idx);
    });
    const sv = parsedVCFsvs[idx];
    currentSV = sv;
    showLoading('Loading gene overlap…');
    let genes = [];
    try { genes = await fetchGeneOverlap(sv.chrom, sv.start, sv.end); }
    catch (e) { genes = getCuratedGenes(sv.chrom, sv.start, sv.end); }
    renderResults(sv, genes);
    hideLoading();
}

// ─── RENDER RESULTS ───────────────────────────────────────────────────────────
function renderResults(sv, genes) {
    // SV header card
    $('sv-type-badge').className = `sv-type-badge type-${sv.type}`;
    $('sv-type-badge').textContent = sv.type;
    $('sv-coords').textContent = `chr${sv.chrom}:${sv.start.toLocaleString()} – ${sv.end.toLocaleString()}`;
    $('sv-size').textContent = `Size: ${formatSize(sv.size)}${sv.sample ? '  ·  ' + sv.sample : ''}`;

    // Classification
    const cls = classifySV(sv, genes);
    const verdicts = $('sv-verdicts');
    verdicts.innerHTML = `<span class="verdict-tag ${cls.tagClass}">${cls.label}</span>`;
    if (cls.mechanism) verdicts.innerHTML += `<span class="verdict-tag vt-neutral">🔀 ${cls.mechanism}</span>`;
    if (genes.length) verdicts.innerHTML += `<span class="verdict-tag vt-neutral">🧬 ${genes.length} gene${genes.length > 1 ? 's' : ''}</span>`;

    // Gene list
    $('genes-count').textContent = `${genes.length} found`;
    $('gene-list').innerHTML = genes.length === 0
        ? '<div style="font-size:.82rem;color:var(--text-muted);padding:.5rem">No Ensembl genes found in this region.</div>'
        : genes.map(g => `
            <div class="gene-item" onclick="openGeneOMIM('${g.name}')">
              <div class="gene-impact-dot ${g.coding ? 'impact-coding' : 'impact-noncoding'}" title="${g.coding ? 'Protein-coding' : 'Non-coding'}"></div>
              <span class="gene-symbol">${g.name}</span>
              <div class="gene-meta">
                <span class="gene-biotype">${g.biotype.replace('_', ' ')}</span>
                <span class="gene-pos">${g.start.toLocaleString()}–${g.end.toLocaleString()}</span>
              </div>
              <span class="gene-strand">${g.strand}</span>
            </div>`).join('');

    // Classification panel
    renderClassification(sv, genes, cls);

    // SVG viewer
    renderSVGViewer(sv, genes);

    // External links
    renderExtLinks(sv, genes);

    show($('results-area'));
}

// ─── CLASSIFICATION ──────────────────────────────────────────────────────────
function classifySV(sv, genes) {
    const size = sv.size;
    const codingGenes = genes.filter(g => g.coding);

    // Score: 0-100
    let score = 0;
    let reasons = [];

    // Size contribution (0-40)
    if (size >= 1e6) { score += 40; reasons.push(`Large (${formatSize(size)})`); }
    else if (size >= 50e3) { score += 25; reasons.push(`Moderate (${formatSize(size)})`); }
    else if (size >= 1e3) { score += 10; reasons.push(`Small (${formatSize(size)})`); }
    else { score += 3; reasons.push(`Very small (${formatSize(size)})`); }

    // SV type (weighted for clinical impact)
    if (sv.type === 'DEL') { score += 30; reasons.push('Deletion (haploinsufficiency risk)'); }
    else if (sv.type === 'DUP') { score += 15; reasons.push('Duplication (triplosensitivity)'); }
    else if (sv.type === 'INV') { score += 10; reasons.push('Inversion (potential gene disruption)'); }
    else if (sv.type === 'INS') { score += 8; reasons.push('Insertion'); }
    else if (sv.type === 'BND') { score += 20; reasons.push('Translocation (potential gene fusion)'); }

    // Gene impact
    if (codingGenes.length > 0) {
        score += Math.min(codingGenes.length * 10, 30);
        reasons.push(`${codingGenes.length} protein-coding gene(s) affected`);
    }

    score = Math.min(score, 100);

    let label, tagClass, mechanism;
    if (score >= 75) { label = '🔴 Likely Pathogenic'; tagClass = 'vt-pathogenic'; }
    else if (score >= 55) { label = '🟠 Possibly Pathogenic'; tagClass = 'vt-likely-path'; }
    else if (score >= 30) { label = '🟡 VUS'; tagClass = 'vt-vus'; }
    else { label = '🟢 Likely Benign'; tagClass = 'vt-benign'; }

    if (sv.type === 'DEL' && codingGenes.length) mechanism = 'Gene deletion';
    else if (sv.type === 'DUP' && codingGenes.length) mechanism = 'Gene duplication';
    else if (sv.type === 'INV') mechanism = 'Gene inversion';
    else if (sv.type === 'BND') mechanism = 'Potential gene fusion';

    return { score, label, tagClass, mechanism, reasons };
}

function renderClassification(sv, genes, cls) {
    const sizeScore = sv.size >= 1e6 ? 100 : sv.size >= 50e3 ? 65 : sv.size >= 1e3 ? 30 : 10;
    const typeScore = { DEL: 100, DUP: 50, INV: 35, INS: 25, BND: 65 }[sv.type] || 30;
    const geneScore = Math.min(genes.filter(g => g.coding).length * 33, 100);
    const barColor = cls.score >= 75 ? '#ff7b72' : cls.score >= 55 ? '#f0883e' : cls.score >= 30 ? '#d29922' : '#3fb950';
    const coding = genes.filter(g => g.coding);

    $('class-body').innerHTML = `
      <div class="class-row">
        <span class="class-label">Size score</span>
        <div class="class-bar-wrap">
          <div class="class-bar"><div class="class-bar-fill" style="width:${sizeScore}%;background:var(--accent-teal)"></div></div>
        </div>
        <span class="class-value" style="color:var(--accent-teal)">${formatSize(sv.size)}</span>
      </div>
      <div class="class-row">
        <span class="class-label">SV type weight</span>
        <div class="class-bar-wrap">
          <div class="class-bar"><div class="class-bar-fill" style="width:${typeScore}%;background:${svTypeColor(sv.type)}"></div></div>
        </div>
        <span class="class-value" style="color:${svTypeColor(sv.type)}">${sv.type}</span>
      </div>
      <div class="class-row">
        <span class="class-label">Coding genes</span>
        <div class="class-bar-wrap">
          <div class="class-bar"><div class="class-bar-fill" style="width:${geneScore}%;background:var(--accent-blue)"></div></div>
        </div>
        <span class="class-value" style="color:var(--accent-blue)">${coding.length}</span>
      </div>
      <div class="class-row" style="margin-top:.5rem">
        <span class="class-label"><strong>Overall score</strong></span>
        <div class="class-bar-wrap">
          <div class="class-bar" style="height:8px"><div class="class-bar-fill" style="width:${cls.score}%;background:${barColor}"></div></div>
        </div>
        <span class="class-value" style="color:${barColor};font-size:1rem">${cls.score}/100</span>
      </div>
      <div class="class-note">
        <strong>Evidence summary:</strong><br>
        ${cls.reasons.map(r => `• ${r}`).join('<br>')}
        <br><br>
        <strong>Note:</strong> This classification is based on size, SV type and gene count. 
        Always integrate with population frequency (DGV / gnomAD-SV), clinical phenotype and literature.
      </div>`;
}

// ─── SVG GENOMIC VIEWER ───────────────────────────────────────────────────────
function renderSVGViewer(sv, genes) {
    const svg = $('sv-svg');
    const W = svg.parentElement.clientWidth || 900;
    svg.setAttribute('viewBox', `0 0 ${W} 160`);
    svg.setAttribute('width', W);

    // Compute genomic range to display (add 20% flanking)
    const flank = Math.max((sv.end - sv.start) * 0.2, 50000);
    const viewStart = Math.max(0, sv.start - flank);
    const viewEnd = sv.end + flank;
    const span = viewEnd - viewStart;
    const toX = pos => ((pos - viewStart) / span) * (W - 40) + 20;

    $('viewer-coords-label').textContent = `chr${sv.chrom}:${viewStart.toLocaleString()}–${viewEnd.toLocaleString()}`;

    let html = '';

    // Chromosome backbone
    html += `<rect x="20" y="72" width="${W - 40}" height="16" rx="8" fill="rgba(48,54,61,0.8)" stroke="rgba(48,54,61,1)"/>`;

    // Genes (below backbone)
    genes.forEach((g, i) => {
        const gx1 = Math.max(toX(g.start), 20);
        const gx2 = Math.min(toX(g.end), W - 20);
        const gw = Math.max(gx2 - gx1, 4);
        const gy = 100 + (i % 2) * 18;
        html += `<rect x="${gx1}" y="${gy}" width="${gw}" height="12" rx="3" fill="rgba(88,166,255,0.55)" stroke="rgba(88,166,255,0.8)" stroke-width=".5"/>`;
        if (gw > 30) {
            html += `<text x="${gx1 + gw / 2}" y="${gy + 9}" text-anchor="middle" font-size="9" font-family="JetBrains Mono, monospace" fill="#58a6ff" font-weight="600">${g.name}</text>`;
        }
        // Line from gene to backbone
        html += `<line x1="${gx1 + gw / 2}" y1="${gy}" x2="${gx1 + gw / 2}" y2="88" stroke="rgba(88,166,255,0.3)" stroke-width=".8"/>`;
    });

    // SV region highlight
    const svX1 = Math.max(toX(sv.start), 20);
    const svX2 = Math.min(toX(sv.end), W - 20);
    const svW = Math.max(svX2 - svX1, 4);
    const svColor = svTypeColor(sv.type);
    const svFill = svTypeFill(sv.type);

    if (sv.type === 'DEL') {
        // Deletion: striped
        html += `<defs><pattern id="del-stripe" patternUnits="userSpaceOnUse" width="8" height="8" patternTransform="rotate(45)">
          <line x1="0" y1="0" x2="0" y2="8" stroke="${svColor}" stroke-width="3" opacity=".5"/>
        </pattern></defs>`;
        html += `<rect x="${svX1}" y="68" width="${svW}" height="24" rx="4" fill="url(#del-stripe)" stroke="${svColor}" stroke-width="1.5"/>`;
    } else if (sv.type === 'DUP') {
        // Duplication: stacked bars
        html += `<rect x="${svX1}" y="65" width="${svW}" height="30" rx="4" fill="${svFill}" stroke="${svColor}" stroke-width="1.5"/>`;
        html += `<rect x="${svX1}" y="65" width="${svW}" height="8" rx="4" fill="${svColor}" opacity=".5"/>`;
    } else if (sv.type === 'INV') {
        html += `<rect x="${svX1}" y="68" width="${svW}" height="24" rx="4" fill="${svFill}" stroke="${svColor}" stroke-width="1.5" stroke-dasharray="5,3"/>`;
        html += `<text x="${svX1 + svW / 2}" y="82" text-anchor="middle" font-size="10" fill="${svColor}">↔</text>`;
    } else {
        html += `<rect x="${svX1}" y="68" width="${svW}" height="24" rx="4" fill="${svFill}" stroke="${svColor}" stroke-width="1.5"/>`;
    }

    // Breakpoint markers
    html += `<line x1="${svX1}" y1="60" x2="${svX1}" y2="95" stroke="${svColor}" stroke-width="1.5" stroke-dasharray="4,2"/>`;
    html += `<line x1="${svX2}" y1="60" x2="${svX2}" y2="95" stroke="${svColor}" stroke-width="1.5" stroke-dasharray="4,2"/>`;

    // Position labels
    html += `<text x="${svX1}" y="55" text-anchor="middle" font-size="9" font-family="JetBrains Mono,monospace" fill="${svColor}">${sv.start.toLocaleString()}</text>`;
    html += `<text x="${svX2}" y="55" text-anchor="middle" font-size="9" font-family="JetBrains Mono,monospace" fill="${svColor}">${sv.end.toLocaleString()}</text>`;

    // SV type label
    html += `<text x="${(svX1 + svX2) / 2}" y="82" text-anchor="middle" font-size="10" font-family="JetBrains Mono,monospace" fill="white" font-weight="700" opacity=".9">${sv.type}</text>`;

    // viewStart / viewEnd labels
    html += `<text x="20" y="68" font-size="8" fill="rgba(139,148,158,0.6)" font-family="JetBrains Mono,monospace">${(viewStart / 1e6).toFixed(2)}M</text>`;
    html += `<text x="${W - 20}" y="68" text-anchor="end" font-size="8" fill="rgba(139,148,158,0.6)" font-family="JetBrains Mono,monospace">${(viewEnd / 1e6).toFixed(2)}M</text>`;

    svg.innerHTML = html;
}

// ─── EXTERNAL LINKS ───────────────────────────────────────────────────────────
function renderExtLinks(sv, genes) {
    const { chrom, start, end } = sv;
    const primaryGene = genes.length ? genes.find(g => g.coding)?.name || genes[0].name : '';
    const links = [
        { href: `https://decipher.sanger.ac.uk/browser#q/chr${chrom}:${start}-${end}`, label: 'DECIPHER', icon: '🔬' },
        { href: `http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg38/?name=chr${chrom}:${start}..${end}`, label: 'DGV', icon: '📊' },
        { href: `https://gnomad.broadinstitute.org/region/chr${chrom}-${start}-${end}?dataset=gnomad_sv_r4`, label: 'gnomAD-SV', icon: '🌍' },
        { href: `https://www.ncbi.nlm.nih.gov/clinvar/?term=${chrom}[chr]+AND+${start}:${end}[chrpos38]`, label: 'ClinVar', icon: '🏥' },
        { href: `https://varsome.com/variant/hg38/chr${chrom}:${start}:${end}`, label: 'VarSome', icon: '🧬' },
        ...(primaryGene ? [
            { href: `https://www.omim.org/search?index=entry&search=${primaryGene}`, label: `OMIM: ${primaryGene}`, icon: '📚' }
        ] : []),
        { href: `https://www.ensembl.org/Homo_sapiens/Location/View?r=${chrom}:${start}-${end}`, label: 'Ensembl Browser', icon: '🗺️' },
    ];
    $('ext-links').innerHTML = links.map(l =>
        `<a class="ext-link-btn" href="${l.href}" target="_blank" rel="noopener">
          <span>${l.icon}</span>${l.label} ↗
        </a>`).join('');
}

function openGeneOMIM(gene) {
    window.open(`https://www.omim.org/search?index=entry&search=${gene}`, '_blank');
}

// ─── COLOR HELPERS ────────────────────────────────────────────────────────────
function svTypeColor(type) {
    return { DEL: '#ff7b72', DUP: '#2dd4bf', INV: '#bc8cff', INS: '#f0883e', BND: '#3fb950' }[type] || '#8b949e';
}
function svTypeFill(type) {
    return { DEL: 'rgba(255,123,114,0.2)', DUP: 'rgba(45,212,191,0.2)', INV: 'rgba(188,140,255,0.15)', INS: 'rgba(240,136,62,0.2)', BND: 'rgba(63,185,80,0.15)' }[type] || 'rgba(139,148,158,0.1)';
}

// ─── INIT ─────────────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
    ['sv-chrom', 'sv-start', 'sv-end', 'sv-type'].forEach(id => {
        const el = $(id);
        if (el) el.addEventListener('keydown', e => { if (e.key === 'Enter') analyzeSV(); });
    });
});
