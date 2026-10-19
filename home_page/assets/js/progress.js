// Draws the spin-1 record-over-time chart from the rows of #progress-table: one step line per
// threshold length L0 (the best bound so far for even L >= L0), on a logarithmic axis, with a
// direct label at the end of each line. The table stays the source of truth.
(function () {
  'use strict';
  const fig = document.getElementById('progress-chart');
  const rows = [...document.querySelectorAll('#progress-table tbody tr')];
  if (!fig || !rows.length) return;
  const pts = rows.map(r => ({date: r.dataset.date, value: parseFloat(r.dataset.value), L0: r.dataset.l0}));

  const W = 640, H = 250, left = 58, right = 92, top = 14, bottom = 34;
  const day = 864e5, t = d => Date.parse(d + 'T00:00:00Z');
  const t0 = Math.min(...pts.map(p => t(p.date))) - day, t1 = Math.max(...pts.map(p => t(p.date))) + day;
  const lo = Math.floor(Math.log10(Math.min(...pts.map(p => p.value))));
  const hi = Math.ceil(Math.log10(Math.max(...pts.map(p => p.value))));
  const X = d => left + (W - left - right) * (t(d) - t0) / (t1 - t0);
  const Y = v => top + (H - top - bottom) * (hi - Math.log10(v)) / (hi - lo);

  const NS = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs, text) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text != null) e.textContent = text;
    return e;
  };
  const svg = el('svg', {viewBox: `0 0 ${W} ${H}`, class: 'chart'});

  for (let k = lo; k <= hi; k++) {
    svg.append(el('line', {x1: left, x2: W - right, y1: Y(10 ** k), y2: Y(10 ** k), class: 'grid'}));
    svg.append(el('text', {x: left - 8, y: Y(10 ** k) + 4, 'text-anchor': 'end', class: 'tick'},
      k >= 0 ? String(10 ** k) : (10 ** k).toFixed(-k)));
  }
  const fmt = d => new Date(t(d)).toLocaleDateString('en-GB', {day: 'numeric', month: 'short', timeZone: 'UTC'});
  for (const d of [...new Set(pts.map(p => p.date))])
    svg.append(el('text', {x: X(d), y: H - bottom + 18, 'text-anchor': 'middle', class: 'tick'}, fmt(d)));

  // Longest series first, so the main line is drawn on top.
  const series = [...new Set(pts.map(p => p.L0))].map(L0 => pts.filter(p => p.L0 === L0)
    .sort((a, b) => a.date.localeCompare(b.date) || a.value - b.value));
  series.sort((a, b) => a.length - b.length);
  const labels = [];
  series.forEach((s, i) => {
    const main = i === series.length - 1;
    let best = 0, path = '';
    for (const p of s) if (p.value > best) {
      path += path ? ` H${X(p.date)} V${Y(p.value)}` : `M${X(p.date)} ${Y(p.value)}`;
      best = p.value;
    }
    path += ` H${W - right}`;
    svg.append(el('path', {d: path, class: main ? 'step main' : 'step'}));
    for (const p of s) svg.append(el('circle', {cx: X(p.date), cy: Y(p.value), r: 3, class: main ? 'pt main' : 'pt'}));
    labels.push({y: Y(best), text: 'L ≥ ' + s[0].L0, main});
  });
  // Keep end labels at least 14 px apart.
  labels.sort((a, b) => a.y - b.y);
  for (let i = 1; i < labels.length; i++) labels[i].y = Math.max(labels[i].y, labels[i - 1].y + 14);
  for (const l of labels)
    svg.append(el('text', {x: W - right + 6, y: l.y + 4, class: l.main ? 'end main' : 'end'}, l.text));

  fig.append(svg);
})();
