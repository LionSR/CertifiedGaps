// The interactive views of leaderboard.html, rendered from the one JSON blob the page embeds
// (#lb-data: tiers, families, records, challenges and the Lean records). Without JavaScript the
// page shows its static table instead.
(function () {
  'use strict';
  const D = JSON.parse(document.getElementById('lb-data').textContent);
  const REPO = 'https://github.com/' + D.repo;
  const RANKED = ['lean-unconditional', 'lean-conditional', 'certified'];
  const tier = Object.fromEntries(D.tiers.map(t => [t.key, t]));
  tier['lean-unavailable'] = {key: 'lean-unavailable', glyph: '?', label: 'Lean status not available in this build', meaning: ''};
  const family = Object.fromEntries(D.families.map(f => [f.key, f]));
  const leanByName = Object.fromEntries((D.lean || []).map(r => [r.name, r]));

  // Resolve every record once: the tier of its headline bound, its Lean theorem, its share.
  const leanTier = r => r.lean ? (leanByName[r.lean.name] ? 'lean-' + leanByName[r.lean.name].kind : 'lean-unavailable') : null;
  const records = D.records.map(r => {
    const value = r.bound ? parseFloat(r.bound.decimal) : null;
    const gap = r.numerical_gap ? parseFloat(r.numerical_gap) : null;
    return Object.assign({}, r, {
      fam: family[r.family], value, gap,
      share: value != null && gap ? value / gap : null,
      eff: r.tier.startsWith('lean-') ? leanTier(r) : r.tier,
      leanEff: leanTier(r),
    });
  });
  const ranked = records.filter(r => r.tier !== 'target');
  const targets = records.filter(r => r.tier === 'target');

  // ---------------------------------------------------------------- small helpers
  function h(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v == null || v === false) continue;
      if (k === 'html') e.innerHTML = v; else if (k.startsWith('on')) e.addEventListener(k.slice(2), v); else e.setAttribute(k, v);
    }
    for (const k of kids.flat(Infinity)) if (k != null) e.append(k instanceof Node ? k : String(k));
    return e;
  }
  // MathJax typesetting must not overlap, so the calls are queued.
  let queue = Promise.resolve();
  const typeset = el => { const M = window.MathJax;
    if (M && M.startup && M.startup.promise) queue = queue.then(() => M.startup.promise).then(() => M.typesetPromise([el])).catch(() => {}); };
  const pct = x => (100 * x).toFixed(1) + '%';
  // LaTeX to plain text, for sorting, search and issue titles.
  const SYMBOLS = {delta: 'δ', theta: 'θ', beta: 'β', perp: '⊥', in: '∈', approx: '≈', le: '≤', ge: '≥', sigma: 'σ'};
  const plain = s => String(s).replace(/\$/g, '').replace(/\\([a-zA-Z]+)\s?/g, (m, c) => SYMBOLS[c] ? SYMBOLS[c] + m.slice(c.length + 1) : '')
    .replace(/\\[ ,]/g, ' ').replace(/[{}_]/g, '');
  const prLink = n => n ? h('a', {href: REPO + '/pull/' + n}, '#' + n) : null;
  const linkPRs = s => String(s).replace(/(^|[\s(])#(\d+)\b/g, (m, a, n) => `${a}<a href="${REPO}/pull/${n}">#${n}</a>`);
  const repoLink = p => /^https?:/.test(p) ? p : REPO + '/blob/' + D.branch + '/' + p;
  const badge = (key, short) => h('span', {class: 'tier tier-' + key, title: tier[key].meaning || tier[key].label},
    tier[key].glyph + ' ' + (short ? '' : tier[key].label));
  const issue = (title, id) => REPO + '/issues/new?template=challenge.yml&title=' + encodeURIComponent(title) +
    '&challenge=' + encodeURIComponent(id);
  const lengthsText = r => (r.lengths.parity === 'odd' ? 'odd ' : r.lengths.parity === 'all' ? '' : 'even ') +
    'L ≥ ' + r.lengths.L0 + (r.lengths.boundary === 'open' ? ', open chain' : '');
  const rowLink = r => h('a', {href: '#' + r.id, onclick: e => { e.preventDefault(); reveal(r.id); }},
    r.fam.name + ', ', h('span', {html: r.parameter}));

  // ---------------------------------------------------------------- highlights
  function highlights() {
    const box = document.getElementById('highlights');
    const closest = ranked.filter(r => r.share != null).sort((a, b) => b.share - a.share).slice(0, 5);
    const newest = [...ranked].sort((a, b) => b.date.localeCompare(a.date) || b.value - a.value).slice(0, 5);
    const lean = ranked.filter(r => r.leanEff && r.leanEff !== 'lean-unavailable');
    const count = k => lean.filter(r => r.leanEff === k).length;
    box.append(
      h('section', {class: 'panel'}, h('h3', {}, 'Closest to the numerical gap'),
        h('ol', {}, closest.map(r => h('li', {}, rowLink(r), ': ', h('strong', {}, pct(r.share)))))),
      h('section', {class: 'panel'}, h('h3', {}, 'Latest records'),
        h('ol', {}, newest.map(r => h('li', {}, rowLink(r), ' ', h('span', {class: 'muted'}, r.date, ' ', prLink(r.pr)))))),
      h('section', {class: 'panel'}, h('h3', {}, 'Formalized in Lean 4'),
        h('p', {}, badge('lean-unconditional', true), ' ', count('lean-unconditional'), ' unconditional (seed + recurrence in Lean kernel)'),
        h('p', {}, badge('lean-conditional', true), ' ', count('lean-conditional'), ' conditional on computer-certified seed inequalities'),
        h('p', {class: 'muted'}, `${lean.length} of ${ranked.length} certified records include a Lean 4 theorem.`)));
    typeset(box);
  }

  // ---------------------------------------------------------------- the leaderboard table
  const state = {tiers: new Set(), fams: new Set(), lean: false, q: '', sort: 'share', dir: -1, open: new Set()};
  const COLS = [
    ['model', 'Model and parameters', 'Hamiltonian family and coupling parameters', r => r.fam.name + ' ' + plain(r.parameter)],
    ['value', 'Certified gap ≥', 'Rigorous lower bound on Δ_L (rounded down) across all covered ring lengths', r => r.value],
    ['share', 'Fraction of numerical gap', 'Certified lower bound divided by the numerical benchmark gap', r => r.share],
    ['tier', 'Verification', 'Verification tier of the displayed bound (and companion Lean 4 theorem where present)', r => RANKED.indexOf(r.eff)],
    ['lengths', 'Lengths covered', 'Ring lengths L covered by the theorem', r => r.lengths.L0],
    ['n0', 'Seed (n₀, β₀)', 'Seed ring length n₀ and inverse temperature β₀ of the finite certificate', r => r.seed.n0],
    ['date', 'Date', 'Date the current bound was established', r => r.date],
  ];

  function controls() {
    const box = document.getElementById('board-controls');
    const chip = (label, title, on, toggle) => h('button', {type: 'button', class: 'chip' + (on() ? ' on' : ''),
      'aria-pressed': on() ? 'true' : 'false', title, onclick: e => { toggle(); e.currentTarget.classList.toggle('on', on());
        e.currentTarget.setAttribute('aria-pressed', on() ? 'true' : 'false'); board(); }}, label);
    const toggleIn = (set, k) => () => set.has(k) ? set.delete(k) : set.add(k);
    const fams = D.families.filter(f => ranked.some(r => r.family === f.key));
    box.append(
      h('input', {type: 'search', placeholder: 'Search models and parameters', 'aria-label': 'Search',
        oninput: e => { state.q = e.target.value.toLowerCase(); board(); }}),
      h('div', {class: 'chips', role: 'group', 'aria-label': 'Verification tier'},
        RANKED.map(k => chip(tier[k].glyph + ' ' + tier[k].label, tier[k].meaning, () => state.tiers.has(k), toggleIn(state.tiers, k))),
        chip('Has a Lean 4 theorem', 'Records accompanied by a Lean 4 theorem (51 of 52 records)',
          () => state.lean, () => { state.lean = !state.lean; })),
      h('div', {class: 'chips', role: 'group', 'aria-label': 'Model family'},
        fams.map(f => chip(f.name, f.name, () => state.fams.has(f.key), toggleIn(state.fams, f.key)))));
  }

  function visible() {
    return ranked.filter(r => (!state.tiers.size || state.tiers.has(r.eff) || (r.leanEff && state.tiers.has(r.leanEff))) &&
      (!state.fams.size || state.fams.has(r.family)) && (!state.lean || r.lean) &&
      (!state.q || (r.fam.name + ' ' + plain(r.parameter)).toLowerCase().includes(state.q)));
  }

  function sorted(rows) {
    const key = COLS.find(c => c[0] === state.sort)[3];
    return rows.sort((a, b) => {
      const x = key(a), y = key(b);
      if (x == null || y == null) return x == null ? (y == null ? 0 : 1) : -1;   // missing values last
      return (x < y ? -1 : x > y ? 1 : 0) * state.dir;
    });
  }

  function details(r) {
    const exact = r.bound.r ? `$\\Delta_L > ${r.bound.display || '\\log(' + r.bound.r + ')/(' + r.bound.beta + ')'}$`
      : r.bound.value.length <= 24 ? `$\\Delta_L ${r.bound.strict ? '>' : '\\ge'} ${r.bound.value}$` : `$\\Delta_L \\ge$ ${r.bound.decimal}…`;
    const title = 'Improve: ' + r.fam.name + ', ' + plain(r.parameter);
    return h('div', {class: 'detail'},
      h('p', {}, h('strong', {}, 'Bound. '), h('span', {html: exact}), ' for every ', lengthsText(r), ' (', r.fam.cells, ').'),
      r.lean ? h('p', {}, h('strong', {}, 'Lean 4 theorem. '), badge(r.leanEff), ' ',
        h('a', {href: D.docs + '?pattern=' + encodeURIComponent(r.lean.name) + '#doc'}, h('code', {}, r.lean.name.replace('Leaderboard.Records.', ''))),
        ' proves ', h('span', {html: `$\\Delta_L > \\log(${r.lean.r})/(${r.lean.beta})$`}), ' ≈ ', r.lean.decimal, ' for every even ', lengthsText(r).replace('even ', ''), '.') : null,
      h('p', {}, h('strong', {}, 'Seed. '), 'n₀ = ', r.seed.n0, ' ', r.fam.cells, ', β₀ = ', r.seed.beta0,
        r.seed.text ? h('span', {html: '; ' + r.seed.text}) : null, '.'),
      r.gap ? h('p', {}, h('strong', {}, 'Numerical benchmark. '), r.numerical_gap, ' (', h('span', {html: r.gap_source || r.fam.gap_source}), ').') : null,
      r.notes ? h('p', {html: linkPRs(r.notes)}) : null,
      h('p', {}, h('strong', {}, 'Certificate. '), h('a', {href: repoLink(r.certificate)}, 'certificate data'), ' · ',
        h('a', {href: repoLink(r.check || 'models/' + (r.slug || r.fam.slug) + '/check.jl')}, 'check.jl'),
        ' · established ', r.date, r.pr ? [' in ', prLink(r.pr)] : null, ' by ',
        r.contributors.map((c, i) => [i ? ', ' : '', person(c)])),
      h('p', {}, h('a', {class: 'button', href: issue(title, r.id)}, 'Improve this record')));
  }

  function grid(fam, rows) {
    const rv = [...new Set(rows.map(r => r.params[fam.grid.rows]))], cv = [...new Set(rows.map(r => r.params[fam.grid.cols]))];
    const at = (a, b) => rows.find(r => r.params[fam.grid.rows] === a && r.params[fam.grid.cols] === b);
    return h('div', {class: 'table-wrap'}, h('table', {class: 'grid'},
      h('caption', {}, 'Certified gap lower bound at each ', fam.grid.rows, ' (rows) and ', fam.grid.cols, ' (columns); select a cell for details.'),
      h('thead', {}, h('tr', {}, h('th', {}, fam.grid.rows + ' \\ ' + fam.grid.cols), cv.map(c => h('th', {}, c)))),
      h('tbody', {}, rv.map(a => h('tr', {}, h('th', {}, a), cv.map(b => {
        const r = at(a, b);
        return h('td', {class: 'num'}, r ? h('a', {href: '#' + r.id, onclick: e => { e.preventDefault(); reveal(r.id); }},
          r.value.toFixed(4), ' ', h('span', {class: 'glyph glyph-' + r.eff, title: tier[r.eff].label}, tier[r.eff].glyph)) : '—');
      }))))));
  }

  function board() {
    const rows = sorted(visible());
    document.getElementById('board-count').textContent = `Showing ${rows.length} of ${ranked.length}`;
    const head = document.getElementById('board-family');
    head.replaceChildren();
    if (state.fams.size === 1) {
      const f = family[[...state.fams][0]];
      head.append(h('p', {class: 'family-head'}, h('strong', {}, f.name, '. '), h('span', {html: f.hamiltonian}),
        '. Ring lengths count ', f.cells, '. Numerical benchmark: ', h('span', {html: f.gap_source}), '.'));
      if (f.grid) head.append(grid(f, ranked.filter(r => r.family === f.key)));
    }
    const tbody = h('tbody');
    for (const r of rows) {
      const open = state.open.has(r.id);
      const companionLean = r.leanEff && r.leanEff !== r.eff ? h('span', {
        class: 'tier tier-' + r.leanEff,
        title: `${tier[r.leanEff].label}: proves Δ_L > ${r.lean.decimal}`
      }, tier[r.leanEff].glyph + ' Lean 4') : null;
      tbody.append(h('tr', {id: r.id, class: open ? 'open' : null},
        h('td', {class: 'model'}, h('span', {class: 'muted'}, r.fam.name), h('br'), h('span', {html: r.parameter})),
        h('td', {class: 'num'}, h('strong', {}, r.bound.decimal)),
        h('td', {class: 'share'}, r.share == null ? h('span', {class: 'muted'}, 'no estimate') :
          h('span', {class: 'bar', title: `${r.bound.decimal} of ${r.numerical_gap}`},
            h('span', {class: 'fill', style: `width:${Math.min(100, 100 * r.share).toFixed(1)}%`}), h('span', {class: 'label'}, pct(r.share)))),
        h('td', {}, badge(r.eff), companionLean ? [' ', companionLean] : null),
        h('td', {}, lengthsText(r)),
        h('td', {class: 'num'}, `(${r.seed.n0}, ${r.seed.beta0})`),
        h('td', {class: 'num'}, r.date),
        h('td', {}, h('button', {type: 'button', class: 'more', 'aria-expanded': open ? 'true' : 'false',
          onclick: () => { open ? state.open.delete(r.id) : state.open.add(r.id); board(); }}, open ? 'less' : 'details'))));
      if (open) tbody.append(h('tr', {class: 'detail-row'}, h('td', {colspan: COLS.length + 1}, details(r))));
    }
    const thead = h('thead', {}, h('tr', {}, COLS.map(([k, label, title]) => h('th', {title,
      'aria-sort': state.sort === k ? (state.dir > 0 ? 'ascending' : 'descending') : null},
      h('button', {type: 'button', class: 'sort', onclick: () => {
        state.dir = state.sort === k ? -state.dir : (k === 'model' || k === 'lengths' ? 1 : -1); state.sort = k; board();
      }}, label, state.sort === k ? (state.dir > 0 ? ' ▲' : ' ▼') : ''))), h('th', {}, '')));
    const table = document.getElementById('board-table');
    table.replaceChildren(thead, tbody);
    typeset(head); typeset(table);
  }

  // Open a row: clear filters that hide it, expand it, scroll to it.
  function reveal(id) {
    const r = ranked.find(x => x.id === id);
    if (!r) return;
    if (!visible().includes(r)) { state.tiers.clear(); state.lean = false; state.q = '';
      state.fams = new Set([r.family]);
      document.querySelectorAll('#board-controls .chip').forEach(c => c.classList.toggle('on', c.textContent === r.fam.name));
      document.querySelector('#board-controls input').value = ''; }
    state.open.add(id);
    board();
    history.replaceState(null, '', '#' + id);
    document.getElementById(id).scrollIntoView({block: 'center'});
  }

  // ---------------------------------------------------------------- challenges
  // Cards for the next records: the hand-written list, the open targets, and challenges derived
  // from the rows (in each family, the point whose bound is the smallest share of its numerical
  // gap; points without a Lean theorem; Lean theorems that still assume their starting inequalities).
  function challenges() {
    const groups = [];
    const group = (title, cards) => { if (cards.length) groups.push({title, cards}); };
    const card = (id, title, need, link) => ({id, title, need, link});
    const perFamily = (pred, make) => D.families.map(f => [f, ranked.filter(r => r.family === f.key && pred(r))])
      .filter(([, rows]) => rows.length).map(([f, rows]) => make(f, rows)).filter(Boolean);
    const some = (n, one, many) => n === 1 ? one : many;
    group('Verification and reproduction', D.challenges.map(c => card(c.id, c.title, c.need, c.link)));
    group('Certify open targets', targets.map(r => {
      const e = r.estimate, from = e.n0 ? ' at seed (n₀, β₀) = (' + e.n0 + (e.beta0 ? ', ' + e.beta0 : '') + ')' : '';
      return card(r.id, 'Certify gap: ' + r.fam.name + ', ' + plain(r.parameter), (e.text ? e.text + '. ' : '') +
        (e.float_bound ? 'Floating-point survey estimate: Δ ≥ ' + e.float_bound + from + '.' : (e.text ? '' : 'No seed estimate yet.')));
    }));
    group('Tighten certified bounds', perFamily(r => r.share != null, (f, rows) => {
      if (rows.length < 2) return null;
      const w = rows.reduce((a, b) => (b.share < a.share ? b : a));
      return card(w.id + '-raise', 'Tighten bound: ' + f.name + ', ' + plain(w.parameter),
        `Current certified bound: ${w.bound.decimal} (${pct(w.share)} of the numerical gap ${w.numerical_gap}, the lowest fraction in this family). ` +
        'A larger seed length n₀ or a sharper spectral enclosure would tighten it.');
    }));
    group('Formalize in Lean 4', perFamily(r => !r.lean, (f, rows) => card(f.key + '-lean', 'Formalize in Lean 4: ' + f.name,
      `No Lean 4 theorem yet for ${some(rows.length, 'parameter point ', rows.length + ' parameter points: ')}` +
      rows.map(r => plain(r.parameter)).join('; ') + '. Formalize a seeded Lean 4 theorem from the certified seed inequalities.')));
    group('Prove seed inequalities in Lean 4', perFamily(r => r.leanEff === 'lean-conditional', (f, rows) =>
      card(f.key + '-unconditional', 'Prove seed in Lean 4: ' + f.name,
        `${rows.length} ${some(rows.length, 'Lean 4 theorem in this family is', 'Lean 4 theorems in this family are')} conditional on computer-certified seed purity inequalities. ` +
        'Proving the finite seed inequalities inside the Lean kernel (as done for the Ising chain at g = 1/4, h = 1/20) ' +
        `makes ${some(rows.length, 'it', 'them')} unconditional.`)));
    const box = document.getElementById('challenge-cards');
    for (const g of groups) box.append(h('h3', {}, g.title, ' ', h('span', {class: 'muted'}, '(' + g.cards.length + ')')),
      h('div', {class: 'cards'}, g.cards.map(c => h('article', {class: 'card'},
        h('h4', {}, c.title), h('p', {html: c.need}),
        h('p', {}, h('a', {class: 'button', href: issue(c.title, c.id)}, 'Take this on'),
          c.link ? [' ', h('a', {href: repoLink(c.link)}, 'background')] : null)))));
    document.getElementById('challenge-count').textContent = groups.reduce((n, g) => n + g.cards.length, 0);
    typeset(box);
  }

  // ---------------------------------------------------------------- contributors
  function person(c) {
    return /^[A-Za-z0-9-]+$/.test(c) ? h('a', {href: 'https://github.com/' + c}, '@' + c) : c;
  }
  function contributors() {
    const tally = {};
    for (const r of ranked) for (const c of r.contributors || []) {
      const t = tally[c] || (tally[c] = {records: 0, lean: 0, improved: 0});
      t.records++; if (r.lean) t.lean++; if (r.history && r.history.length) t.improved++;
    }
    const rows = Object.entries(tally).sort((a, b) => b[1].records - a[1].records || a[0].localeCompare(b[0]));
    document.getElementById('contributor-table').append(
      h('thead', {}, h('tr', {}, h('th', {}, 'Contributor'),
        h('th', {title: 'Records where this contributor established the current bound'}, 'Current records'),
        h('th', {title: 'Current records that include a Lean 4 theorem'}, 'With Lean 4 theorem'),
        h('th', {title: 'Current records that improved a previous bound'}, 'Improved previous bound'))),
      h('tbody', {}, rows.map(([c, t]) => h('tr', {}, h('td', {}, person(c)),
        h('td', {class: 'num'}, t.records), h('td', {class: 'num'}, t.lean), h('td', {class: 'num'}, t.improved)))));
  }

  // ---------------------------------------------------------------- timeline
  function timelineGroups() {
    const groups = [];
    for (const f of D.families) {
      const rows = ranked.filter(r => r.family === f.key);
      const point = (r, x) => ({date: x.date, value: parseFloat(x.decimal || x.bound.decimal), tier: x.tier === r.tier ? r.eff : x.tier,
        lengths: x.lengths || lengthsText(r), what: x.method || 'current bound'});
      const series = r => [point(r, r)].concat((r.history || []).map(x => point(r, x)));
      if (rows.length && rows.every(r => !r.params || !Object.keys(r.params).length)) {
        const pts = rows.flatMap(series);
        if (pts.length > 1) groups.push({label: f.name, pts});
      } else rows.filter(r => r.history && r.history.length).forEach(r =>
        groups.push({label: f.name + ', ' + plain(r.parameter), pts: series(r)}));
    }
    return groups;
  }
  function timeline() {
    const groups = timelineGroups();
    if (!groups.length) return;
    const pick = h('select', {'aria-label': 'Model', onchange: e => draw(groups[e.target.value])},
      groups.map((g, i) => h('option', {value: i}, g.label)));
    document.getElementById('timeline-pick').append(h('label', {}, 'Model: ', pick));
    draw(groups[0]);
  }
  function draw(g) {
    const NS = 'http://www.w3.org/2000/svg', W = 640, H = 260, L = 70, R = 90, T = 16, B = 44;
    const s = (tag, attrs, text) => { const e = document.createElementNS(NS, tag);
      for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (text != null) e.textContent = text; return e; };
    const pts = [...g.pts].sort((a, b) => a.date.localeCompare(b.date) || a.value - b.value);
    const t = d => Date.parse(d + 'T00:00:00Z'), day = 864e5;
    const t0 = t(pts[0].date) - 2 * day, t1 = t(pts[pts.length - 1].date) + 2 * day;
    const lo = Math.floor(Math.log10(Math.min(...pts.map(p => p.value)))), hi = Math.ceil(Math.log10(Math.max(...pts.map(p => p.value))));
    const X = d => L + (W - L - R) * (t(d) - t0) / (t1 - t0), Y = v => T + (H - T - B) * (hi - Math.log10(v)) / Math.max(1, hi - lo);
    const svg = s('svg', {viewBox: `0 0 ${W} ${H}`, role: 'img', class: 'timeline',
      'aria-label': `Best certified gap over time for ${g.label}; the table below lists the points.`});
    for (let k = lo; k <= hi; k++) {
      svg.append(s('line', {x1: L, x2: W - R, y1: Y(10 ** k), y2: Y(10 ** k), class: 'grid'}));
      svg.append(s('text', {x: L - 6, y: Y(10 ** k) + 4, 'text-anchor': 'end', class: 'tick'}, k >= 0 ? String(10 ** k) : (10 ** k).toFixed(-k)));
    }
    const dates = [...new Set(pts.map(p => p.date))];
    dates.forEach(d => svg.append(s('text', {x: X(d), y: H - B + 16, 'text-anchor': 'middle', class: 'tick'}, d.slice(5))));
    svg.append(s('text', {x: (L + W - R) / 2, y: H - 6, 'text-anchor': 'middle', class: 'axis'}, 'Date (2026)'));
    svg.append(s('text', {x: 14, y: (T + H - B) / 2, transform: `rotate(-90 14 ${(T + H - B) / 2})`, 'text-anchor': 'middle', class: 'axis'},
      'Certified gap ≥ (log scale)'));
    // The step line follows the best bound so far.
    let best = 0, path = '';
    for (const p of pts) if (p.value > best) {
      path += path ? ` H${X(p.date)} V${Y(p.value)}` : `M${X(p.date)} ${Y(p.value)}`; best = p.value;
    }
    path += ` H${W - R}`;
    svg.append(s('path', {d: path, class: 'step'}));
    for (const p of pts) {
      const c = s('circle', {cx: X(p.date), cy: Y(p.value), r: 5, class: 'pt glyph-' + p.tier});
      c.append(s('title', {}, `${p.date}: ${p.value} (${p.lengths}; ${tier[p.tier] ? tier[p.tier].label : p.tier}; ${p.what})`));
      svg.append(c);
    }
    svg.append(s('text', {x: W - R + 6, y: Y(best) + 4, class: 'endlabel'}, String(best)));
    document.getElementById('timeline-chart').replaceChildren(svg);
    document.getElementById('timeline-table').replaceChildren(
      h('thead', {}, h('tr', {}, h('th', {}, 'Date'), h('th', {}, 'Certified gap ≥'), h('th', {}, 'Lengths covered'), h('th', {}, 'Method'))),
      h('tbody', {}, pts.slice().reverse().map(p => h('tr', {}, h('td', {}, p.date), h('td', {class: 'num'}, p.value),
        h('td', {}, p.lengths), h('td', {}, p.what)))));
  }

  // ---------------------------------------------------------------- start
  highlights();
  controls();
  board();
  challenges();
  contributors();
  timeline();
  // A link may name a row (#<id>) or a model family (#family-<key>).
  const hash = decodeURIComponent(location.hash.slice(1));
  if (ranked.some(r => r.id === hash)) reveal(hash);
  else if (hash.startsWith('family-') && family[hash.slice(7)]) {
    state.fams = new Set([hash.slice(7)]);
    document.querySelectorAll('#board-controls .chip').forEach(c => c.classList.toggle('on', c.textContent === family[hash.slice(7)].name));
    board();
    document.getElementById('board').scrollIntoView();
  }
})();
