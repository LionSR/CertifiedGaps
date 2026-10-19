// Draws the spin-1 record-over-time chart from the rows of #progress-table. The bounds are
// placed in the order they appeared (equal steps, dates below), and a staircase follows the
// best bound so far. Every bound in the table holds for all even L >= 2304, so the staircase is
// the best known bound for large rings. A hollow point is a bound that did not raise it.
(function () {
  'use strict';
  const fig = document.getElementById('progress-chart');
  const rows = [...document.querySelectorAll('#progress-table tbody tr')];
  if (!fig || !rows.length) return;
  const pts = rows.map(r => ({date: r.dataset.date, text: r.dataset.value,
    value: parseFloat(r.dataset.value), L0: r.dataset.l0}))
    .sort((a, b) => a.date.localeCompare(b.date) || a.value - b.value);

  const W = 640, H = 240, left = 52, right = 24, top = 22, bottom = 40;
  const step = (W - left - right) / pts.length;
  const X = i => left + step * (i + 0.5);
  const vmax = Math.ceil(Math.max(...pts.map(p => p.value)) / 0.02) * 0.02;
  const Y = v => top + (H - top - bottom) * (1 - v / vmax);

  // Labels cut the stored decimal (already rounded down) to four significant digits; rounding
  // could raise a lower bound.
  const sig4 = s => s.replace(/^(0\.0*[1-9]\d{0,3}|[1-9]\d*\.\d{0,3}).*$/, '$1');
  const NS = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs, text) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text != null) e.textContent = text;
    return e;
  };
  const svg = el('svg', {viewBox: `0 0 ${W} ${H}`, class: 'chart'});

  for (let v = 0; v <= vmax + 1e-9; v += 0.02) {
    svg.append(el('line', {x1: left, x2: W - right, y1: Y(v), y2: Y(v), class: 'grid'}));
    svg.append(el('text', {x: left - 8, y: Y(v) + 4, 'text-anchor': 'end', class: 'tick'}, v.toFixed(2)));
  }
  // One date label under each run of bounds from the same day.
  const fmt = d => new Date(Date.parse(d + 'T00:00:00Z'))
    .toLocaleDateString('en-GB', {day: 'numeric', month: 'short', timeZone: 'UTC'});
  for (let i = 0; i < pts.length;) {
    let j = i;
    while (j + 1 < pts.length && pts[j + 1].date === pts[i].date) j++;
    svg.append(el('text', {x: (X(i) + X(j)) / 2, y: H - bottom + 20, 'text-anchor': 'middle', class: 'tick'}, fmt(pts[i].date)));
    i = j + 1;
  }

  let best = 0, path = '';
  const marks = pts.map((p, i) => {
    const up = p.value > best;
    if (up) { path += path ? ` H${X(i)} V${Y(p.value)}` : `M${X(i)} ${Y(p.value)}`; best = p.value; }
    return {p, i, up};
  });
  path += ` H${W - right}`;
  svg.append(el('path', {d: path, class: 'step main'}));
  for (const {p, i, up} of marks) {
    svg.append(el('circle', {cx: X(i), cy: Y(p.value), r: 3.5, class: up ? 'pt main filled' : 'pt'}));
    svg.append(el('text', {x: X(i), y: Y(p.value) - 9, 'text-anchor': 'middle', class: up ? 'end main' : 'end'},
      sig4(p.text)));
  }
  fig.append(svg);
})();
