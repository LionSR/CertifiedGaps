// Animated chapters of the spacetime-bootstrap explainer. Each chapter draws one frame at time t
// on a 960x540 logical canvas, so the page players and the rendered film share every picture.
(function (root) {
  'use strict';
  const W = 960, H = 540;
  const clamp = (x, a = 0, b = 1) => Math.max(a, Math.min(b, x));
  const seg = (t, a, b) => clamp((t - a) / (b - a));
  const ease = x => { x = clamp(x); return x < .5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2; };
  const lerp = (a, b, s) => a + (b - a) * s;
  // Tufte: one serif face throughout, thin strokes, no rounded panels, direct labels.
  const SERIF = '"EB Garamond", Palatino, "DejaVu Serif", serif';

  // ---------- drawing kit ----------
  function kit(ctx, th, texImg) {
    const K = {ctx, th, W, H, seg, ease, lerp, clamp};
    K.text = (s, x, y, o = {}) => {
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.fillStyle = o.color || th.ink;
      ctx.font = `${o.italic ? 'italic ' : ''}${o.weight || 400} ${Math.round((o.size || 16) * 1.12)}px ${SERIF}`; ctx.textAlign = o.align || 'left';
      ctx.textBaseline = o.base || 'alphabetic'; ctx.fillText(s, x, y); ctx.restore();
    };
    K.tex = (src, x, y, o = {}) => {
      const im = texImg(src, o.color || th.ink, o.size || 20); if (!im) return;
      const xx = o.align === 'center' ? x - im.w / 2 : o.align === 'right' ? x - im.w : x;
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.drawImage(im.img, xx, y - im.h / 2, im.w, im.h); ctx.restore();
    };
    K.line = (x1, y1, x2, y2, o = {}) => {
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.strokeStyle = o.color || th.ink; ctx.lineWidth = Math.max(.8, (o.w || 2) * .6);
      if (o.dash) ctx.setLineDash(o.dash); ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); ctx.restore();
    };
    K.arrow = (x1, y1, x2, y2, o = {}) => {
      const a = Math.atan2(y2 - y1, x2 - x1), hd = (o.head || 10) * .75;
      K.line(x1, y1, x2 - Math.cos(a) * hd * .6, y2 - Math.sin(a) * hd * .6, o);
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.fillStyle = o.color || th.ink; ctx.beginPath();
      ctx.moveTo(x2, y2); ctx.lineTo(x2 - hd * Math.cos(a - .42), y2 - hd * Math.sin(a - .42));
      ctx.lineTo(x2 - hd * Math.cos(a + .42), y2 - hd * Math.sin(a + .42)); ctx.closePath(); ctx.fill(); ctx.restore();
    };
    K.rect = (x, y, w, h, o = {}) => {
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1;
      if (o.fill) { ctx.fillStyle = o.fill; ctx.fillRect(x, y, w, h); }
      if (o.stroke) { ctx.strokeStyle = o.stroke; ctx.lineWidth = Math.max(.8, (o.w || 2) * .5); if (o.dash) ctx.setLineDash(o.dash); ctx.strokeRect(x, y, w, h); }
      ctx.restore();
    };
    K.dot = (x, y, r, color, alpha = 1) => { ctx.save(); ctx.globalAlpha *= alpha; ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, r * .8, 0, 7); ctx.fill(); ctx.restore(); };
    K.circle = (x, y, r, o = {}) => { ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.strokeStyle = o.color || th.ink; ctx.lineWidth = Math.max(.8, (o.w || 2) * .6); if (o.dash) ctx.setLineDash(o.dash); ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.stroke(); ctx.restore(); };
    K.path = (pts, o = {}) => {
      if (pts.length < 2) return; ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.strokeStyle = o.color || th.ink; ctx.lineWidth = Math.max(1, (o.w || 2.5) * .65);
      ctx.lineJoin = 'round'; if (o.dash) ctx.setLineDash(o.dash); ctx.beginPath(); ctx.moveTo(pts[0][0], pts[0][1]);
      for (const p of pts.slice(1)) ctx.lineTo(p[0], p[1]); ctx.stroke(); ctx.restore();
    };
    // bars: values in [0,1], baseline at y+h
    K.bars = (x, y, w, h, vals, color, o = {}) => {
      const n = vals.length, gap = o.gap ?? 8, bw = (w - gap * (n - 1)) / n;
      K.line(x - 4, y + h, x + w + 4, y + h, {color: th.mute, w: 1.5, alpha: o.alpha});
      vals.forEach((v, i) => {
        K.rect(x + i * (bw + gap) + bw * .2, y + h - v * h, bw * .6, v * h, {fill: color, alpha: (o.alpha ?? 1) * (i === 0 || o.flat ? .9 : .45)});
        if (o.labels && v > .03) K.text(v.toFixed(2), x + i * (bw + gap) + bw / 2, y + h - v * h - 7, {size: 12, align: 'center', color: th.ink, alpha: o.alpha});
      });
    };
    // flat spacetime torus: cols x rows cells of size c; dashed identified edges
    K.torus = (x, y, cols, rows, c, o = {}) => {
      const w = cols * c, h = rows * c, al = o.alpha ?? 1;
      K.rect(x, y, w, h, {fill: o.fill || th.surface, alpha: al});
      for (let i = 1; i < cols; i++) K.line(x + i * c, y, x + i * c, y + h, {color: th.line, w: 1, alpha: al});
      for (let j = 1; j < rows; j++) K.line(x, y + j * c, x + w, y + j * c, {color: th.line, w: 1, alpha: al});
      K.line(x, y, x, y + h, {color: th.space, dash: [5, 4], alpha: al}); K.line(x + w, y, x + w, y + h, {color: th.space, dash: [5, 4], alpha: al});
      K.line(x, y, x + w, y, {color: th.time, dash: [5, 4], alpha: al}); K.line(x, y + h, x + w, y + h, {color: th.time, dash: [5, 4], alpha: al});
    };
    K.axes = (x, y, w, h, o = {}) => { K.line(x, y + h, x + w, y + h, {color: th.mute, w: 1.5, alpha: o.alpha}); K.line(x, y, x, y + h, {color: th.mute, w: 1.5, alpha: o.alpha}); };
    // a tensor drawn as a labeled box centered at (x, y); its legs are drawn separately as lines
    K.box = (x, y, w, h, label, o = {}) => {
      const col = o.color || th.proved;
      K.rect(x - w / 2, y - h / 2, w, h, {fill: th.surface, alpha: o.alpha});
      K.rect(x - w / 2, y - h / 2, w, h, {fill: o.fill || th.provedF, stroke: col, w: 2, alpha: o.alpha});
      if (label) K.text(label, x, y + 6, {align: 'center', color: col, size: o.size || 16, alpha: o.alpha});
    };
    return K;
  }

  const CH = [];
  const chapter = c => CH.push(c);

  // 1 ------------------------------------------------------------------
  chapter({
    id: 'gap', title: 'What is a gap?', dur: 15,
    beats: [[0, 'A chain of L spins has a ladder of energy levels. The gap Δ is the first step up.'],
            [4.5, 'Haldane (1983): for spin 1 the gap stays open as the chain grows; for spin ½ it closes.'],
            [10, 'Goal: one number c > 0 with Δ(L) ≥ c for every ring length L.']],
    draw(K, t) {
      const {th} = K, a = seg(t, 0, .8);
      const lx = 70, base = 400, sc = 150, E = [0, .9, .95, 1, 1.3, 1.45, 1.6, 1.8];
      E.forEach((e, i) => K.line(lx, base - e * sc, lx + 170, base - e * sc, {color: i ? th.mute : th.ink, w: i ? 2 : 4, alpha: a * seg(t, .2 + i * .1, .6 + i * .1)}));
      const g = seg(t, 1.5, 2.3);
      K.arrow(lx + 200, base, lx + 200, base - .9 * sc + 2, {color: th.proved, w: 2.5, alpha: g}); K.arrow(lx + 200, base - .9 * sc, lx + 200, base - 2, {color: th.proved, w: 2.5, alpha: g});
      K.text('Δ', lx + 214, base - .45 * sc + 8, {size: 26, color: th.proved, alpha: g});
      K.text('E₀', lx - 10, base + 6, {align: 'right', color: th.mute, alpha: a}); K.text('energy', lx, 110, {color: th.mute, alpha: a, size: 14});
      // plot of Δ_L versus L
      const px = 380, py = 110, pw = 520, ph = 290, b = seg(t, 4.5, 5.3);
      const X = L => px + L / 64 * pw, Y = d => py + ph - d / 1.6 * ph;
      K.axes(px, py, pw, ph, {alpha: b});
      K.text('ring length L', px + pw, py + ph + 28, {align: 'right', color: th.mute, alpha: b, size: 14});
      K.text('gap Δ(L)', px - 10, py - 8, {color: th.mute, alpha: b, size: 14});
      const shown = clamp((t - 5) / 4) * 31;
      for (let i = 0; i <= 30; i++) {
        if (i > shown) break; const L = 4 + 2 * i;
        K.dot(X(L), Y(5.5 / L), 4, th.time, b); K.dot(X(L), Y(.41 + 9 / (L * L)), 4, th.space, b);
      }
      K.text('spin ½: closes', X(40), Y(.32), {color: th.time, alpha: seg(t, 7, 8)});
      K.text('spin 1: stays near 0.41', X(30), Y(.6), {color: th.space, alpha: seg(t, 7.5, 8.5)});
      K.text('schematic', px + pw, py - 8, {size: 12, color: th.mute, alpha: b, align: 'right'});
      const c = seg(t, 10, 11);
      K.rect(px, Y(.25), pw, Y(0) - Y(.25), {fill: th.provedF, alpha: c});
      K.line(px, Y(.25), px + pw, Y(.25), {color: th.proved, dash: [7, 5], w: 2.5, alpha: c});
      K.text('c', px + pw + 10, Y(.25) + 6, {size: 22, color: th.proved, alpha: c});
      K.text('every spin-1 dot must stay above the dashed line', px, py + ph + 56, {color: th.proved, alpha: seg(t, 11.5, 12.5), size: 15});
    }
  });

  // 2 ------------------------------------------------------------------
  chapter({
    id: 'rate', title: 'A gap is a decay rate', dur: 14,
    beats: [[0, 'At inverse temperature β, level m has probability e^(−βEₘ)/Z.'],
            [4, 'Cooling down empties the excited levels. The ratio of the first two weights is exactly e^(−βΔ).'],
            [9.5, 'So it is enough to show: excited weight ≤ e^(−cβ) on every ring. Then Δ ≥ c.']],
    draw(K, t) {
      const {th} = K, E = [0, 1, 1, 1, 1.6, 2, 2.4], beta = lerp(.15, 4.5, ease(seg(t, 2, 9)));
      const w = E.map(e => Math.exp(-beta * e)), Z = w.reduce((x, y) => x + y);
      K.bars(70, 110, 330, 280, w.map(v => v / Z), th.time, {labels: true, alpha: seg(t, 0, .8)});
      K.text(`β = ${beta.toFixed(2)}`, 70, 440, {size: 18, alpha: seg(t, 0, .8)});
      K.text('levels  E = 0, 1, 1, 1, 1.6, 2, 2.4', 70, 470, {color: th.mute, size: 14, alpha: seg(t, 0, .8)});
      const px = 500, py = 110, pw = 400, ph = 280, b = seg(t, 4, 4.8), Y = v => py - v / 4.5 * ph;
      K.axes(px, py, pw, ph, {alpha: b});
      const pts = []; for (let bb = 0; bb <= beta; bb += .05) pts.push([px + bb / 4.5 * pw, py + bb / 4.5 * ph]);
      K.path(pts, {color: th.proved, w: 3, alpha: b});
      K.text('log (w₁ / w₀)', px + 8, py - 12, {color: th.mute, size: 14, alpha: b});
      K.text('β →', px + pw, py + ph + 26, {color: th.mute, align: 'right', size: 14, alpha: b});
      K.text('slope = −Δ', px + pw * .55, py + ph * .45, {color: th.proved, size: 18, alpha: seg(t, 6, 7)});
      K.tex('\\frac{w_1}{w_0}=e^{-\\beta\\Delta}', 700, 470, {align: 'center', size: 26, alpha: seg(t, 9.5, 10.5)});
    }
  });

  // A ------------------------------------------------------------------
  chapter({
    id: 'imag', title: 'Imaginary time', dur: 22,
    beats: [[0, 'A quantum state is a list of numbers cₘ, one per energy level. Draw each as an arrow: its length is |cₘ|, its angle is the phase. The bars below are the probabilities |cₘ|².'],
            [5, 'Real time, e^(−iHt): arrow m turns at speed Eₘ. Lengths never change, so the probabilities never change. Nothing settles.'],
            [9.5, 'Now replace t by −iτ. The same formula gives e^(−τH): no turning at all. Arrow m shrinks by the factor e^(−τEₘ); high energies shrink fastest.'],
            [14, 'Rescale so the probabilities add to 1. The ground state takes over, and the first excited level loses like e^(−τΔ). The gap is the speed.'],
            [18.5, 'At τ = β the shrink factors are exactly the Boltzmann weights e^(−βE). Low temperature and long imaginary time are the same thing.']],
    draw(K, t) {
      const {th} = K, E = [0, 1, 1.5, 2.2], c = [.45, .6, .45, .5], cols = [th.space, th.time, th.proved, th.ok];
      const real = t >= 5 && t < 9.5, tr = real ? (t - 5) * 1.3 : 0, tau = t >= 9.5 ? lerp(0, 3, ease(seg(t, 10, 13.5))) : 0;
      const renorm = seg(t, 14, 15.5), R = 50, y = 160;
      const amp = c.map((v, m) => v * Math.exp(-tau * E[m])), n2 = amp.reduce((s, v) => s + v * v, 0), n20 = c.reduce((s, v) => s + v * v, 0);
      const scale = lerp(1 / Math.sqrt(n20), 1 / Math.sqrt(n2), renorm);
      K.text(real ? 'real time: arrows turn, lengths fixed' : t >= 9.5 ? `imaginary time τ = ${tau.toFixed(2)}${renorm > .5 ? ', rescaled' : ''}` : 'a state, level by level', 70, 60, {size: 18});
      E.forEach((e, m) => {
        const x = 130 + m * 140, len = amp[m] * scale * R * 1.3, ang = -e * tr;
        K.circle(x, y, R, {color: th.line, w: 1.5});
        K.arrow(x, y, x + len * Math.cos(ang), y - len * Math.sin(ang), {color: cols[m], w: 3.5, head: 9});
        K.text(m ? `level ${m}` : 'ground', x, y - R - 12, {align: 'center', color: cols[m], size: 14});
        K.text(`E = ${e}`, x, y + R + 22, {align: 'center', color: th.mute, size: 13});
        const p = Math.pow(amp[m] * scale, 2), bh = 110;
        K.rect(x - 22, 400 - p * bh, 44, p * bh, {fill: cols[m], alpha: .85});
        K.text(p.toFixed(2), x, 392 - p * bh, {align: 'center', color: th.mute, size: 12});
        if (t >= 18.5) K.text(`e^(−βE) = ${Math.exp(-3 * e).toFixed(3)}`, x, 432, {align: 'center', color: th.proved, size: 13, alpha: seg(t, 18.8, 19.5)});
      });
      K.line(90, 400, 610, 400, {color: th.line, w: 1.5}); K.text('prob.', 60, 404, {color: th.mute, size: 13, align: 'right'});
      // ratio panel
      const b = seg(t, 14, 14.8);
      if (b > 0) {
        const gx = 690, gy = 110, gw = 220, gh = 200, Xp = v => gx + v / 3 * gw, Yp = v => gy - v / 3 * gh * 1;
        K.axes(gx, gy, gw, gh, {alpha: b});
        const pts = []; for (let v = 0; v <= tau + 1e-9; v += .05) pts.push([Xp(v), gy + v / 3 * gh]);
        K.path(pts, {color: th.time, w: 3, alpha: b});
        K.text('log (level 1 / ground)', gx, gy - 12, {color: th.mute, size: 13, alpha: b});
        K.text('slope −Δ', gx + gw * .55, gy + gh * .45, {color: th.time, size: 15, alpha: b});
        K.text('τ →', gx + gw, gy + gh + 20, {align: 'right', color: th.mute, size: 13, alpha: b});
      }
      K.tex('e^{-\\tau H}\\psi=\\sum_m c_m\\,e^{-\\tau E_m}|m\\rangle', 790, 400, {align: 'center', size: 20, alpha: seg(t, 10.5, 11.5)});
    }
  });

  // B ------------------------------------------------------------------
  // Matrix products as sums over routes, shown before the paths of A^M: the reason weights multiply.
  function routes(K, t) {
    const {th} = K, xs = [140, 320, 500], ys = [180, 310], g = seg(t, 6.5, 7.3);
    [['start i', 0], ['middle j', 1], ['end k', 2]].forEach(([s, i]) => K.text(s, xs[i], 130, {align: 'center', color: th.mute, size: 14}));
    [[0, th.space, ['A₁₁', '2'], ['B₁₁', '2']], [1, th.time, ['A₁₂', '1'], ['B₂₁', '1']]].forEach(([j, col, l1, l2], q) => {
      const al = seg(t, 1.5 + q * .8, 2.2 + q * .8), my = (ys[0] + ys[j]) / 2 - 12;
      K.line(xs[0], ys[0], xs[1], ys[j], {color: col, w: 3, alpha: al}); K.line(xs[1], ys[j], xs[2], ys[0], {color: col, w: 3, alpha: al});
      [[l1, (xs[0] + xs[1]) / 2], [l2, (xs[1] + xs[2]) / 2]].forEach(([l, x]) => {
        K.text(l[0], x, my, {align: 'center', color: col, size: 15, alpha: al * (1 - g)}); K.text(l[1], x, my, {align: 'center', color: col, size: 17, alpha: al * g}); });
    });
    xs.forEach((x, c) => ys.forEach((y, s) => { K.dot(x, y, 6, th.ink);
      if (c === 0) K.text(String(s + 1), x - 14, y + 5, {align: 'right', color: th.mute, size: 14});
      else K.text(String(s + 1), x, s ? y + 24 : y - 24, {align: 'center', color: th.mute, size: 14}); }));
    K.tex('(AB)_{ik}=\\sum_j A_{ij}\\,B_{jk}', 760, 140, {align: 'center', size: 24, alpha: seg(t, .5, 1.3)});
    K.text('one product per route, added over all routes', 610, 195, {color: th.mute, size: 15, alpha: seg(t, 3, 3.8)});
    // the same rule as a tensor diagram: a joined leg is a summed index
    const dg = seg(t, 4, 4.8) * (1 - g), y = 270;
    K.line(620, y, 820, y, {color: th.ink, w: 2, alpha: dg});
    K.box(670, y, 40, 30, 'A', {alpha: dg}); K.box(770, y, 40, 30, 'B', {alpha: dg});
    [['i', 628], ['j', 720], ['k', 812]].forEach(([s, x]) => K.text(s, x, y - 10, {align: 'center', color: th.mute, size: 14, italic: true, alpha: dg}));
    K.text('joined leg j = sum over j', 620, 320, {color: th.proved, size: 15, alpha: dg});
    K.tex('A=B=\\begin{pmatrix}2&1\\\\1&3\\end{pmatrix}', 760, 270, {align: 'center', size: 22, alpha: g});
    K.text('route via j = 1:   2 × 2 = 4', 610, 345, {color: th.space, size: 17, alpha: seg(t, 8, 8.7)});
    K.text('route via j = 2:   1 × 1 = 1', 610, 377, {color: th.time, size: 17, alpha: seg(t, 9, 9.7)});
    K.text('(A²)₁₁ = 4 + 1 = 5', 610, 420, {size: 21, alpha: seg(t, 10, 10.7)});
  }
  function pathsBody(K, t) {
    const {th} = K, A1 = 1 - seg(t, 12.2, 12.8);
    const ROWS = [150, 260], lab = ['↑', '↓'];
    if (A1 > 0) {
      const M = 8, x0 = 110, dx = 58, path = [0, 0, 1, 1, 1, 0, 0, 0, 0];
      for (let s = 0; s < M; s++) for (let a = 0; a < 2; a++) for (let b = 0; b < 2; b++) {
        const flip = a !== b, emph = t >= 6;
        K.line(x0 + s * dx, ROWS[a], x0 + (s + 1) * dx, ROWS[b], {color: flip ? th.proved : th.mute, w: emph ? (flip ? 1 : 2) : 1, dash: emph && flip ? [4, 4] : null, alpha: A1 * .45});
      }
      const prog = clamp((t - 1.5) / 3.5) * M;
      const pts = []; for (let s = 0; s <= Math.min(M, Math.floor(prog) + 1); s++) { const ss = Math.min(s, prog); const s0 = Math.floor(ss), f = ss - s0; const yA = ROWS[path[s0]], yB = ROWS[path[Math.min(M, s0 + 1)]]; pts.push([x0 + ss * dx, lerp(yA, yB, f)]); if (s > prog) break; }
      K.path(pts, {color: th.space, w: 4, alpha: A1});
      for (let s = 0; s <= M; s++) for (let a = 0; a < 2; a++) K.dot(x0 + s * dx, ROWS[a], 5, th.ink, A1);
      lab.forEach((l, a) => K.text(l, x0 - 26, ROWS[a] + 6, {align: 'center', size: 20, alpha: A1}));
      K.text('step →', x0 + M * dx, ROWS[1] + 40, {align: 'right', color: th.mute, size: 13, alpha: A1});
      const nflip = path.slice(0, M).reduce((s, v, i) => s + (v !== path[i + 1] ? 1 : 0), 0);
      K.text(t < 6 ? 'weight of the red path = product of the entries it uses' : `red path: ${M - nflip} stays × 1, ${nflip} flips × dτ  →  weight dτ${nflip === 2 ? '²' : ''}`, x0, 340, {size: 16, alpha: A1 * seg(t, 2, 3)});
      K.text('(Aᴹ)ᵢⱼ = sum of the weights of all paths from i to j', x0, 372, {size: 16, color: th.proved, alpha: A1 * seg(t, 4, 5)});
      K.tex('A=\\begin{pmatrix}1&d\\tau\\\\ d\\tau&1\\end{pmatrix}', 780, 200, {align: 'center', size: 24, alpha: A1 * seg(t, 6.5, 7.5)});
      K.text('stay: weight 1', 700, 290, {color: th.ink, size: 15, alpha: A1 * seg(t, 7, 8)});
      K.text('flip: weight dτ', 700, 314, {color: th.proved, size: 15, alpha: A1 * seg(t, 7.5, 8.5)});
      K.text('e^(−βH) = Aᴹ,  M = β / dτ', 700, 350, {color: th.ink, size: 15, alpha: A1 * seg(t, 9, 10)});
    }
    const B = seg(t, 12.8, 13.4) * (1 - seg(t, 22.7, 23.3));
    if (B > 0) {
      const x0 = 110, x1 = 850, Xt = f => lerp(x0, x1, f);
      const hist = [[], [.22, .6], [.15, .4, .7, .85]];
      hist.forEach((fl, h) => {
        const yU = 120 + h * 120, yD = yU + 46, al = B * seg(t, 13.4 + h * .8, 14 + h * .8);
        let up = true; const pts = [[Xt(0), yU]];
        for (const f of fl) { pts.push([Xt(f), up ? yU : yD]); up = !up; pts.push([Xt(f), up ? yU : yD]); }
        pts.push([Xt(1), up ? yU : yD]);
        K.line(x0, yU, x1, yU, {color: th.line, w: 1, alpha: al}); K.line(x0, yD, x1, yD, {color: th.line, w: 1, alpha: al});
        K.path(pts, {color: th.space, w: 4, alpha: al});
        fl.forEach(f => K.dot(Xt(f), (yU + yD) / 2, 5, th.proved, al));
        K.text('↑', x0 - 22, yU + 6, {align: 'center', alpha: al}); K.text('↓', x0 - 22, yD + 6, {align: 'center', alpha: al});
        K.text(`k = ${fl.length}:  weight βᵏ/k! = ${(Math.pow(2, fl.length) / [1, 1, 2, 6, 24][fl.length]).toFixed(2)} (β = 2)`, x1, yD + 26, {align: 'right', color: th.mute, size: 13, alpha: al});
      });
      K.text('imaginary time τ: 0 → β', x0, 92, {color: th.time, size: 14, alpha: B});
      const c = seg(t, 19, 20);
      K.line(x0, 100, x0, 470, {color: th.time, dash: [6, 4], w: 2, alpha: B * c}); K.line(x1, 100, x1, 470, {color: th.time, dash: [6, 4], w: 2, alpha: B * c});
      K.text('τ = 0 and τ = β are the same moment: start state = end state', W / 2, 500, {align: 'center', color: th.time, size: 16, alpha: B * c});
    }
    const C = seg(t, 23.3, 24);
    if (C > 0) {
      const beta = 2, terms = [0, 2, 4, 6, 8].map(k => 2 * Math.pow(beta, k) / [1, 1, 2, 6, 24, 120, 720, 5040, 40320][k]);
      const gx = 160, gy = 110, gw = 420, gh = 300, mx = 8, X = i => gx + 30 + i * 78, Y = v => gy + gh - v / mx * gh;
      K.axes(gx, gy, gw, gh, {alpha: C});
      let run = 0; const nshow = Math.floor(lerp(1, 5.99, seg(t, 24, 26.5)));
      terms.forEach((v, i) => { if (i >= nshow) return; K.rect(X(i), Y(run + v), 50, Y(run) - Y(run + v), {fill: th.space, alpha: C * .85}); K.text(`k=${2 * i}`, X(i) + 25, gy + gh + 20, {align: 'center', color: th.mute, size: 12, alpha: C}); run += v; });
      K.line(gx, Y(2 * Math.cosh(beta)), gx + gw, Y(2 * Math.cosh(beta)), {color: th.proved, dash: [6, 4], w: 2, alpha: C});
      K.text(`2 cosh 2 = ${(2 * Math.cosh(2)).toFixed(3)}`, gx + gw, Y(2 * Math.cosh(beta)) - 8, {align: 'right', color: th.proved, size: 14, alpha: C});
      K.text(`running sum = ${run.toFixed(3)}`, gx, gy - 14, {size: 15, alpha: C});
      K.tex('Z=\\sum_{k\\ \\mathrm{even}}\\frac{2\\beta^k}{k!}=2\\cosh\\beta', 780, 260, {align: 'center', size: 24, alpha: C * seg(t, 25, 26)});
    }
  }
  chapter({
    id: 'paths', title: 'A sum over histories', dur: 41,
    beats: [[0, 'Start with how matrices multiply: (AB)ᵢₖ = Σⱼ Aᵢⱼ Bⱼₖ. To go from i to k, pass through a middle state j, multiply the two entries on the way, and add over every j.'],
            [6.5, 'With numbers: for A = [[2, 1], [1, 3]], (A²)₁₁ = A₁₁A₁₁ + A₁₂A₂₁ = 2·2 + 1·1 = 5, one product per route 1 → j → 1.'],
            [13, 'Repeat for M factors: an entry of Aᴹ is a sum over paths through M steps. Each path picks one entry of A per step, and its weight is the product of those entries.'],
            [19, 'One spin that can flip, H = −σˣ. For a tiny step dτ, A = e^(−dτ H) ≈ [[1, dτ], [dτ, 1]]: staying costs 1, flipping costs dτ. And e^(−βH) = Aᴹ with M = β/dτ.'],
            [25.5, 'Make the steps tiny. A path becomes a world line that mostly stays and flips a few times. With k flips its weight is dτᵏ, and there are about (β/dτ)ᵏ/k! places for them: βᵏ/k! in total.'],
            [32, 'The trace keeps only paths that end where they started. So k is even, and the time axis closes into a circle.'],
            [36, 'Add all paths: Σ 2βᵏ/k! over even k = 2 cosh β, the exact answer. This is the Euclidean path integral: quantum mechanics as a sum over histories.']],
    draw(K, t) {
      const a = 1 - seg(t, 12.4, 13), b = seg(t, 13, 13.6);
      if (a > 0) { K.ctx.save(); K.ctx.globalAlpha *= a; routes(K, t); K.ctx.restore(); }
      if (b > 0) { K.ctx.save(); K.ctx.globalAlpha *= b; pathsBody(K, t - 13); K.ctx.restore(); }
    }
  });

  // 3 ------------------------------------------------------------------
  chapter({
    id: 'torus', title: 'The partition function is a torus', dur: 15,
    beats: [[0, 'Z = Tr e^(−βH). Draw the chain across and imaginary time up: a sheet.'],
            [4, 'The trace makes the last time slice equal the first: imaginary time is a circle.'],
            [8.5, 'The chain is a ring, so space is a circle too. Every quantity in the proof lives on this torus.']],
    draw(K, t) {
      const {th, ctx} = K;
      const m = ease(seg(t, 4, 7.5)) + ease(seg(t, 8.5, 12));
      const yaw = -0.55 + (t > 12 ? (t - 12) * .35 : 0) + .08 * Math.sin(t * .6);
      STB.drawTorus(ctx, {cx: W / 2, cy: 285, scale: lerp(250, 440, seg(m, 1, 2)), m, yaw, pitch: lerp(.35, .9, seg(m, 0, 2)),
        cols: 16, rows: 8, colHi: -1, nu: 64, nv: 32, lineW: 1.5, shadeAmt: .35,
        colors: {space: th.space, time: th.time, fill: th.torusFill, hi: th.proved}});
      const a = 1 - seg(t, 3.5, 4.2);
      K.text('space: L sites →', W / 2 - 170, 470, {color: th.space, size: 18, alpha: a});
      K.text('imaginary time β ↑', W / 2 + 250, 300, {color: th.time, size: 18, alpha: a});
      K.text('time closes', W / 2, 500, {align: 'center', color: th.time, size: 18, alpha: seg(t, 5, 5.6) * (1 - seg(t, 8, 8.5))});
      K.text('space closes', W / 2, 500, {align: 'center', color: th.space, size: 18, alpha: seg(t, 9.5, 10.1)});
    }
  });

  // C ------------------------------------------------------------------
  // Classical warm-up: the Ising transfer matrix as a tensor network (boxes are matrices, legs are
  // indices, a joined leg is a summed index), contracted one spin at a time down to Tr T^L.
  const ISP = [1, 1, -1, 1, -1, 1], SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹', SUB = '₀₁₂₃₄₅₆₇₈₉';
  chapter({
    id: 'ising', title: 'Transfer matrices: the Ising ring', dur: 38,
    beats: [[0, 'Warm-up: a classical Ising ring. Each site carries a spin s = ±1, and each pair of neighbors contributes a weight e^(K s s′): e^K if they agree, e^(−K) if not.'],
            [7, 'A configuration’s weight is the product of its bond weights, and Z adds up the weights of all 2ᴸ configurations.'],
            [12.5, 'Store the bond weight as a 2×2 matrix T(s, s′) = e^(K s s′), drawn as a box with two legs, one per index. The spins become the legs between boxes.'],
            [19.5, 'A leg joining two boxes means: sum over that index. Σ_{s₂} T(s₁, s₂) T(s₂, s₃) = (T²)(s₁, s₃) is matrix multiplication. This is how T passes the spin on to the next site.'],
            [28, 'Join every leg: the row of boxes becomes Tᴸ. The last leg wraps around to the first, and joining those two is a trace: Z = Tr Tᴸ.']],
    draw(K, t) {
      const {th} = K, n = 6, x0 = 150, dx = 120, y = 190, X = i => x0 + i * dx;
      const bl = i => X(i) + 36, br = i => i < n - 1 ? X(i + 1) - 36 : X(i) + 84, bc = i => (bl(i) + br(i)) / 2;
      const tm = seg(t, 12.5, 13.5), A0 = seg(t, .2, .8);
      const k = seg(t, 22, 23.2) + seg(t, 28.5, 29.3) + seg(t, 29.5, 30.3) + seg(t, 30.5, 31.3) + seg(t, 31.5, 32.3);
      const full = Math.min(5, Math.floor(k)), fr = full < 5 ? k - full : 0;
      const alive = m => m === 0 ? 1 : 1 - clamp(k - (m - 1));
      // the wrap-around bond and the legs
      const wr = seg(t, 32.6, 33.4);
      K.path([[br(n - 1), y], [X(n - 1) + 110, y], [X(n - 1) + 110, y + 60], [X(0) - 50, y + 60], [X(0) - 50, y], [X(0) - 26 * (1 - tm), y]],
        {color: wr > 0 ? th.space : th.proved, w: 2 + wr, dash: [6, 4], alpha: A0 * seg(t, 1.5, 2.2)});
      K.line(X(0), y, br(n - 1), y, {color: th.ink, w: 2, alpha: tm});
      for (let i = 0; i < n; i++) {
        const up = ISP[i] > 0, al = A0 * (1 - tm);
        K.circle(X(i), y, 26, {color: th.ink, w: 2, alpha: al});
        K.arrow(X(i), y + (up ? 14 : -14), X(i), y + (up ? -14 : 14), {color: up ? th.space : th.time, w: 3, alpha: al, head: 8});
        const ag = ISP[i] * ISP[(i + 1) % n] > 0;
        K.text(ag ? 'e^K' : 'e^−K', bc(i), y + 6, {align: 'center', color: ag ? th.ok : th.fail, size: 15, alpha: (1 - tm) * seg(t, 1.5 + i * .3, 2 + i * .3)});
        K.dot(X(i), y, 5, th.ink, tm * alive(i));
        K.text(`s${SUB[i + 1]}`, X(i), y + 30, {align: 'center', color: th.mute, size: 15, alpha: tm * alive(i)});
      }
      // boxes: bonds 0..full merged into one box, the rest drawn singly
      const R = full < 5 ? lerp(br(full), br(full + 1), ease(fr)) : br(5), cnt = Math.min(6, full + 1 + (fr > .5 ? 1 : 0));
      for (let m = full + 1; m < n; m++) K.box(bc(m), y, br(m) - bl(m), 30, 'T', {alpha: tm * (m === full + 1 ? 1 - ease(fr) : 1)});
      K.box((bl(0) + R) / 2, y, R - bl(0), 30, cnt === 1 ? 'T' : cnt === 6 ? 'T⁶ = Tᴸ' : 'T' + SUP[cnt], {alpha: tm});
      // beat 0 and 1: weights and their product
      K.text('agree: e^K        disagree: e^−K', W / 2, 330, {align: 'center', color: th.mute, size: 16, alpha: seg(t, 3, 3.8) * (1 - seg(t, 6.6, 7))});
      const b1 = 1 - seg(t, 12, 12.5);
      K.text('this configuration:  e^K · e^−K · e^−K · e^−K · e^−K · e^K = e^(−2K)', W / 2, 330, {align: 'center', size: 17, alpha: seg(t, 7.3, 8) * b1});
      K.text('Z = the sum of such products over all 2⁶ = 64 configurations', W / 2, 368, {align: 'center', size: 17, color: th.proved, alpha: seg(t, 9, 9.8) * b1});
      // beat 2: the matrix and the product of its entries
      K.text('box = matrix,  leg = index', W / 2, 330, {align: 'center', color: th.mute, size: 16, alpha: seg(t, 15, 15.8) * (1 - seg(t, 19.2, 19.6))});
      K.tex('T=\\begin{pmatrix}e^{K}&e^{-K}\\\\ e^{-K}&e^{K}\\end{pmatrix}', 240, 420, {align: 'center', size: 20, alpha: seg(t, 13.5, 14.3)});
      K.tex('\\text{weight}=T(s_1,s_2)\\,T(s_2,s_3)\\cdots T(s_6,s_1)', 640, 420, {align: 'center', size: 20, alpha: seg(t, 16, 16.8)});
      // beat 3: one joined leg is one matrix product
      K.rect(bl(0) - 8, y - 26, br(1) - bl(0) + 16, 52, {stroke: th.space, w: 2, dash: [5, 4], alpha: seg(t, 19.8, 20.4) * (1 - seg(t, 23.2, 23.8))});
      K.tex('\\sum_{s_2}T(s_1,s_2)\\,T(s_2,s_3)=(T^2)(s_1,s_3)', W / 2, 330, {align: 'center', size: 20, alpha: seg(t, 20.2, 21) * (1 - seg(t, 27.4, 27.9))});
      K.text('joined leg = summed index = matrix multiplication', W / 2, 482, {align: 'center', color: th.proved, size: 16, alpha: seg(t, 23.5, 24.2)});
      // beat 4: the wrap-around leg is a trace
      K.tex('Z=\\sum_{s_1}(T^L)(s_1,s_1)=\\operatorname{Tr}T^L', W / 2, 330, {align: 'center', size: 24, color: th.space, alpha: seg(t, 33, 34)});
    }
  });

  // The quantum chain: events from expanding e^(-beta H), one site's time-ordered product, the
  // per-site tensor W and its trace X, and the same tensor-network picture as for Ising.
  // Site 3's events are the example of 'rpwhere': a = {z at 0.2, x at 0.6}, b = {y at 0.4}.
  const CEV = [[0, [[.15, 'z'], [.8, 'x']]], [1, [[.2, 'z'], [.6, 'x']]], [2, [[.4, 'y']]], [3, [[.25, 'x'], [.55, 'y']]], [4, [[.35, 'y'], [.88, 'z']]], [5, [[.12, 'x'], [.7, 'z']]]];
  // a row of four boxes joined by legs, the last leg wrapping around to the first
  function ringRow(K, x, y, lab, legs, col, al) {
    const {th} = K;
    K.path([[x + 270, y], [x + 284, y], [x + 284, y + 24], [x - 14, y + 24], [x - 14, y], [x, y]], {color: th.ink, w: 2, alpha: al});
    K.line(x, y, x + 270, y, {color: th.ink, w: 2, alpha: al});
    for (let i = 0; i < 4; i++) K.box(x + 38 + i * 64, y, 34, 26, lab, {color: col, alpha: al});
    legs.forEach((s, i) => K.text(s, x + 6 + i * 64, y - 10, {align: 'center', color: th.mute, size: 14, alpha: al}));
  }
  chapter({
    id: 'chain', title: 'Histories of a chain', dur: 61,
    beats: [[0, 'Quantum chain. As for the single spin, cut β into tiny steps dτ: e^(−βH) ≈ (1 − dτ H)ᴹ. In each step either nothing happens (weight 1) or one term of −H acts (weight dτ).'],
            [7.5, 'For Heisenberg, −S_j·S_{j+1} = Σ_α G_α ⊗ G_α, with G_α = iS^α real 3×3 matrices. So “one term acts” means: pick a bond and a type α ∈ {x, y, z}, and apply G_α to both sites of that bond.'],
            [15.5, 'What the G’s are: in the basis |x⟩, |y⟩, |z⟩, G_x, G_y, G_z are real antisymmetric 3×3 matrices. G_z turns |x⟩ into −|y⟩ and |y⟩ into |x⟩ and kills |z⟩: a quarter turn about the z axis. For example G_z G_z = diag(−1, −1, 0), whose trace is −2.'],
            [22.5, 'That is an event: a rung across one bond at one imaginary time, labeled by its type. A history is the set of all events; each event carries weight dτ, exactly like a flip of the single spin.'],
            [30, 'Follow site 3. Its events arrive from its left bond (list a, gold) and its right bond (list b, blue): z at τ = 0.2, y at 0.4, x at 0.6. Reading upward, its matrices multiply in time order, later ones on the left.'],
            [37.5, 'As a tensor, W(a, b) has legs s′ at τ = 0 and s at τ = β: W = ⟨s| G_x G_y G_z |s′⟩. Imaginary time is a circle, so join the two time legs (s = s′, a trace): X(a, b) = Tr(G_x G_y G_z) = 1.'],
            [45.5, 'Other cases: no events gives Tr 1 = 3; one z-event on each side gives Tr G_z G_z = −2. The weight of a whole history is a trace over the chain, and it factorizes into one such trace per site.'],
            [53, 'So a history weighs X(a₆, a₁) X(a₁, a₂) ⋯: boxes X on sites, joined by legs that carry bond histories. It is the Ising diagram with T replaced by X, so Z = Tr Xᴸ. This X is the transfer matrix of the proof.']],
    draw(K, t) {
      // a 7 s pause at 15.5 shows the G matrices; the history itself is frozen meanwhile
      const gm = seg(t, 15.5, 16.1) * (1 - seg(t, 21.9, 22.5));
      chainBody(K, t < 15.5 ? t : t < 22.5 ? 15.5 : t - 7, gm);
      gPanel(K, gm);
    }
  });

  function chainBody(K, t, gm) {
      const {th} = K, N = 6, x0 = 110, dx = 88, y0 = 80, y1 = 400, Hh = y1 - y0, X = i => x0 + i * dx, Yt = f => y1 - f * Hh;
      const B = seg(t, .2, .8), hi = seg(t, 23, 23.8) * (1 - seg(t, 45.4, 46)), P1 = (1 - seg(t, 22.4, 23)) * (1 - gm);
      const P2 = seg(t, 23, 23.5) * (1 - seg(t, 38, 38.5)), P3 = seg(t, 38.5, 39) * (1 - seg(t, 45.4, 46)), P4 = seg(t, 46, 46.6);
      // the time steps of the expansion
      for (let q = 1; q < 16; q++) K.line(X(0) - 20, Yt(q / 16), X(N - 1) + 20, Yt(q / 16), {color: th.line, w: 1, alpha: B * seg(t, 1, 2) * (1 - seg(t, 15, 16))});
      K.rect(X(2) - 44, y0 - 6, 88, Hh + 12, {fill: th.spaceF, alpha: hi});
      for (let i = 0; i < N; i++) {
        const me = i === 2 && hi > 0;
        K.line(X(i), y0, X(i), y1, {color: me ? th.ink : th.mute, w: me ? 3 : 2, alpha: B});
        K.text(`${i + 1}`, X(i), y1 + 22, {align: 'center', color: th.mute, size: 13, alpha: B});
        K.box(X(i), y1 + 52, 28, 22, 'X', {size: 14, alpha: seg(t, 46.5 + i * .2, 47 + i * .2)});
        K.text(`a${SUB[i + 1]}`, X(i) + dx / 2, y0 - 14, {align: 'center', color: th.mute, size: 15, alpha: seg(t, 47.5 + i * .2, 48 + i * .2)});
      }
      K.text('site', X(0) - 30, y1 + 22, {align: 'right', color: th.mute, size: 13, alpha: B});
      K.text('τ ↑', X(0) - 30, y0 + 10, {align: 'right', color: th.time, size: 14, alpha: B});
      CEV.forEach(([bd, evs], k) => evs.forEach(([f, ty], q) => {
        const ex = bd === 3 && q === 0, at = ex ? 11 : 16 + (k * 2 + q) * .25, al = B * seg(t, at, at + .5), yy = Yt(f);
        const col = hi > 0 && bd === 1 ? th.proved : hi > 0 && bd === 2 ? th.time : th.space;
        if (bd < N - 1) K.line(X(bd), yy, X(bd + 1), yy, {color: col, w: 3, alpha: al});
        else { K.line(X(bd), yy, X(bd) + 40, yy, {color: col, w: 3, alpha: al, dash: [5, 3]}); K.line(X(0) - 40, yy, X(0), yy, {color: col, w: 3, alpha: al, dash: [5, 3]}); }
        K.dot(X(bd), yy, 4, col, al); K.dot(X((bd + 1) % N), yy, 4, col, al);
        K.text(ty, bd < N - 1 ? X(bd) + dx / 2 : X(bd) + 22, yy - 7, {align: 'center', color: col, size: 15, italic: true, alpha: al});
        if (ex) { const e = al * P1; K.text('G_x', X(bd) - 7, yy + 5, {align: 'right', size: 14, alpha: e}); K.text('G_x', X(bd + 1) + 7, yy + 5, {size: 14, alpha: e}); }
      }));
      K.text('a', X(1) + dx / 2, y0 - 14, {align: 'center', color: th.proved, size: 22, alpha: hi * (1 - seg(t, 46, 46.5))});
      K.text('b', X(2) + dx / 2, y0 - 14, {align: 'center', color: th.time, size: 22, alpha: hi * (1 - seg(t, 46, 46.5))});
      const lp = seg(t, 34, 34.8) * (1 - seg(t, 45.4, 46));
      K.path([[X(2), y0], [X(2), y0 - 30], [X(2) + 30, y0 - 30], [X(2) + 30, y1 + 34], [X(2), y1 + 34], [X(2), y1]], {color: th.space, w: 2, dash: [4, 3], alpha: lp});
      // panel 1: where events come from
      K.tex('e^{-\\beta H}\\approx(1-d\\tau\\,H)^{M}', 770, 100, {align: 'center', size: 22, alpha: P1 * seg(t, .5, 1.3)});
      K.text('each step: nothing (× 1),', 620, 152, {size: 16, alpha: P1 * seg(t, 2.5, 3.2)});
      K.text('or one term of −H (× dτ)', 620, 176, {size: 16, alpha: P1 * seg(t, 3.5, 4.2)});
      K.tex('-\\mathbf S_j\\cdot\\mathbf S_{j+1}=\\sum_{\\alpha=x,y,z}G_\\alpha\\otimes G_\\alpha', 770, 232, {align: 'center', size: 20, alpha: P1 * seg(t, 8, 8.8)});
      K.text('one term = a bond and a type α:', 620, 288, {size: 16, alpha: P1 * seg(t, 10, 10.7)});
      K.text('G_α on both of its sites', 620, 312, {size: 16, alpha: P1 * seg(t, 10.5, 11.2)});
      K.text('event: a rung at one time, type x, y or z', 620, 362, {size: 16, color: th.space, alpha: P1 * seg(t, 16, 16.7)});
      K.text('history: all events, weight dτ each', 620, 388, {size: 16, color: th.space, alpha: P1 * seg(t, 18, 18.7)});
      // panel 2: site 3's product, the tensor W, and its trace X
      K.text('site 3, from the bottom up:', 610, 100, {color: th.mute, size: 15, alpha: P2 * seg(t, 24, 24.7)});
      [['τ = 0.2', 'G_z', 'from a', th.proved], ['τ = 0.4', 'G_y', 'from b', th.time], ['τ = 0.6', 'G_x', 'from a', th.proved]].forEach(([a, g, s, col], i) => {
        const al = P2 * seg(t, 25 + i * .8, 25.6 + i * .8), yy = 130 + i * 26;
        K.text(a, 620, yy, {color: th.mute, size: 15, alpha: al}); K.text(g, 700, yy, {size: 16, alpha: al}); K.text(s, 750, yy, {color: col, size: 15, alpha: al}); });
      K.text('product, later on the left:  G_x G_y G_z', 610, 220, {size: 16, alpha: P2 * seg(t, 28, 28.7)});
      const wv = P2 * seg(t, 31, 31.8), wx = 690, wy = 320;
      K.line(wx, wy - 50, wx, wy + 50, {color: th.ink, w: 2, alpha: wv});
      K.line(wx - 60, wy, wx, wy, {color: th.proved, w: 2.5, alpha: wv}); K.line(wx, wy, wx + 60, wy, {color: th.time, w: 2.5, alpha: wv});
      K.box(wx, wy, 44, 36, 'W', {alpha: wv});
      K.text('s', wx + 8, wy - 40, {size: 15, italic: true, alpha: wv * (1 - lp)}); K.text('s′', wx + 8, wy + 50, {size: 15, italic: true, alpha: wv * (1 - lp)});
      K.text('a', wx - 64, wy + 5, {align: 'right', color: th.proved, size: 16, alpha: wv}); K.text('b', wx + 64, wy + 5, {color: th.time, size: 16, alpha: wv});
      K.path([[wx, wy - 50], [wx, wy - 60], [wx + 84, wy - 60], [wx + 84, wy + 60], [wx, wy + 60], [wx, wy + 50]], {color: th.space, w: 2, dash: [4, 3], alpha: P2 * lp});
      K.text('= ⟨s| G_x G_y G_z |s′⟩', wx + 100, wy + 5, {size: 16, alpha: wv * seg(t, 32, 32.7) * (1 - lp)});
      K.text('joined: s = s′, a trace', wx + 100, wy + 5, {size: 16, color: th.space, alpha: P2 * lp});
      K.text('X(a, b) = Tr(G_x G_y G_z) = 1', 610, 430, {size: 21, alpha: P2 * seg(t, 35.5, 36.3)});
      // panel 3: more values of X, and the factorization over sites
      const rows = [['a', 'b', 'X(a, b)'], ['none', 'none', 'Tr 1 = 3'], ['z', 'none', 'Tr G_z = 0'], ['z', 'z', 'Tr G_zG_z = −2'], ['z, x', 'y', 'Tr G_xG_yG_z = 1']];
      rows.forEach((r, i) => { const yy = 110 + i * 34, al = P3 * seg(t, 39 + i * .5, 39.6 + i * .5);
        K.text(r[0], 610, yy, {color: i ? th.proved : th.mute, size: 16, alpha: al}); K.text(r[1], 690, yy, {color: i ? th.time : th.mute, size: 16, alpha: al});
        K.text(r[2], 770, yy, {color: i ? th.ink : th.mute, size: 16, alpha: al});
        if (i === 0) K.line(610, yy + 10, 920, yy + 10, {color: th.line, w: 1, alpha: al}); });
      K.text('entries can be negative: that is fine', 610, 285, {color: th.mute, size: 15, alpha: P3 * seg(t, 41.5, 42.2)});
      K.text('trace over the chain = one trace per site:', 610, 335, {size: 16, alpha: P3 * seg(t, 42.5, 43.2)});
      K.tex('\\text{weight}=\\prod_j X(a_{j-1},a_j)', 760, 380, {align: 'center', size: 20, alpha: P3 * seg(t, 43.5, 44.3)});
      // panel 4: the same diagram as Ising
      K.text('Ising: spins on the legs', 610, 110, {color: th.mute, size: 15, alpha: P4});
      ringRow(K, 630, 150, 'T', ['s₁', 's₂', 's₃', 's₄'], th.proved, P4 * seg(t, 46.5, 47.2));
      K.text('quantum chain: bond histories on the legs', 610, 230, {color: th.mute, size: 15, alpha: P4 * seg(t, 48, 48.6)});
      ringRow(K, 630, 270, 'X', ['a₁', 'a₂', 'a₃', 'a₄'], th.space, P4 * seg(t, 48.5, 49.2));
      K.text('each joined leg: sum over its history', 610, 345, {color: th.proved, size: 16, alpha: P4 * seg(t, 50, 50.7)});
      K.text('Z = Tr Xᴸ', 760, 420, {align: 'center', size: 30, color: th.space, alpha: seg(t, 51, 52)});
  }

  // the three spin-1 generators in the Cartesian basis
  function gPanel(K, gm) {
    if (gm <= 0) return;
    const {th} = K, mats = ['\\begin{pmatrix}0&0&0\\\\0&0&1\\\\0&-1&0\\end{pmatrix}', '\\begin{pmatrix}0&0&-1\\\\0&0&0\\\\1&0&0\\end{pmatrix}', '\\begin{pmatrix}0&1&0\\\\-1&0&0\\\\0&0&0\\end{pmatrix}'];
    K.text('G_α = iS^α in the basis |x⟩, |y⟩, |z⟩:', 610, 96, {size: 16, alpha: gm});
    ['G_x', 'G_y', 'G_z'].forEach((g, a) => { K.text(g, 645 + a * 105, 130, {align: 'center', size: 16, alpha: gm}); K.tex(mats[a], 645 + a * 105, 180, {align: 'center', size: 15, alpha: gm}); });
    K.text('G_z: |x⟩ → −|y⟩,  |y⟩ → |x⟩,  |z⟩ → 0', 610, 262, {size: 16, alpha: gm});
    K.text('a quarter turn about the z axis', 610, 288, {color: th.mute, size: 15, alpha: gm});
    K.text('G_z G_z = diag(−1, −1, 0),   Tr = −2', 610, 334, {size: 16, alpha: gm});
    K.text('real, antisymmetric: Gᵀ = −G', 610, 372, {color: th.mute, size: 15, alpha: gm});
  }

  // 4 ------------------------------------------------------------------
  chapter({
    id: 'slice', title: 'Two ways to slice the torus', dur: 16,
    beats: [[0, 'Slice across: each row is a step of time evolution. Summing over energies gives Z = Σ e^(−βE).'],
            [5, 'Slice up: each column is one site carried along the ring by a transfer matrix X. Then Z = Tr Xᴸ = Σ λᴸ.'],
            [10.5, 'Same number, two sums. On a long ring the largest λ wins, just as at low temperature the lowest E wins.']],
    draw(K, t) {
      const {th} = K, x = 70, y = 120, c = 40, cols = 8, rows = 6;
      K.torus(x, y, cols, rows, c);
      K.text('L sites', x + cols * c / 2, y + rows * c + 28, {align: 'center', color: th.space, size: 15});
      K.text('β', x - 18, y + rows * c / 2 + 6, {align: 'center', color: th.time, size: 18});
      if (t < 5) { const r = Math.floor(seg(t, .5, 4.5) * rows * 1.999) % rows; K.rect(x, y + r * c, cols * c, c, {fill: th.timeF, stroke: th.time, w: 2, alpha: seg(t, .3, .8)}); }
      else if (t < 10.5) { const k = Math.floor(seg(t, 5.3, 10) * cols * 1.999) % cols; K.rect(x + k * c, y, c, rows * c, {fill: th.spaceF, stroke: th.space, w: 2, alpha: seg(t, 5, 5.5)}); }
      K.tex('Z=\\sum_m e^{-\\beta E_m}', 640, 160, {align: 'center', size: 26, color: th.time, alpha: seg(t, 1, 2)});
      K.tex('Z=\\operatorname{Tr}X^L=\\sum_i\\lambda_i^{\\,L}', 640, 240, {align: 'center', size: 26, color: th.space, alpha: seg(t, 6, 7)});
      const L = Math.round(lerp(1, 16, ease(seg(t, 11.5, 15.5)))), lam = [1, .75, .6, .5, .4];
      const p = lam.map(l => Math.pow(l, L)), s = p.reduce((u, v) => u + v), b = seg(t, 10.5, 11.3);
      K.bars(470, 300, 340, 150, p.map(v => v / s), th.space, {alpha: b, labels: true});
      K.text(`share of Z from each λ, ring length L = ${L}`, 470, 482, {color: th.mute, size: 14, alpha: b});
    }
  });

  // 5 ------------------------------------------------------------------
  chapter({
    id: 'purity', title: 'Purity, seen as squares', dur: 15,
    beats: [[0, 'A state with probabilities p₀, p₁, … has purity P = Σ p². Draw each p as a square: P is their total area.'],
            [5, 'Every square fits under the height of the largest one, p₀. So P ≤ p₀.'],
            [9.5, 'Purity near 1 therefore forces one state to carry almost all the probability.']],
    draw(K, t) {
      const {th} = K, x0 = 90, base = 450, U = 380;
      const k = ease(seg(t, 9.5, 13.5)), p = [lerp(.45, .9, k), lerp(.3, .05, k), lerp(.15, .03, k), lerp(.1, .02, k)];
      const cols = [th.space, th.time, th.proved, th.ok];
      let x = x0;
      p.forEach((v, i) => { const s = v * U; K.rect(x, base - s, s, s, {fill: cols[i], alpha: .32 * seg(t, .5 + i * .4, 1 + i * .4)}); K.rect(x, base - s, s, s, {stroke: cols[i], w: 2, alpha: seg(t, .5 + i * .4, 1 + i * .4)}); x += s; });
      K.line(x0, base, x0 + U, base, {color: th.ink, w: 2});
      K.text('0', x0, base + 22, {align: 'center', color: th.mute, size: 13}); K.text('1', x0 + U, base + 22, {align: 'center', color: th.mute, size: 13});
      K.text('probabilities laid end to end', x0 + U / 2, base + 22, {align: 'center', color: th.mute, size: 13});
      const r = seg(t, 5, 6);
      K.rect(x0, base - p[0] * U, U, p[0] * U, {stroke: th.ink, w: 2, dash: [7, 5], alpha: r});
      K.text('height p₀', x0 + U + 10, base - p[0] * U + 6, {color: th.ink, alpha: r, size: 15});
      const P = p.reduce((s, v) => s + v * v, 0);
      K.text(`P = Σp² = ${P.toFixed(3)}`, 600, 200, {size: 22, alpha: seg(t, 2, 3)});
      K.text(`p₀ = ${p[0].toFixed(3)}`, 600, 240, {size: 22, alpha: r});
      K.text('squares ⊂ dashed box  ⇒  P ≤ p₀', 600, 300, {size: 18, color: th.proved, alpha: seg(t, 6.5, 7.5)});
    }
  });

  // 6 ------------------------------------------------------------------
  // The two density matrices are what is left after cutting the closed network open: cutting one
  // leg of the ring of X's leaves X^n (sigma), cutting every site at one moment leaves e^(-beta H) (rho).
  chapter({
    id: 'two', title: 'Two purities from one torus', dur: 44,
    beats: [[0, 'The torus as a network: one transfer matrix X per site, joined leg to leg around the ring. Joining all the legs gives Zₙ = Tr Xⁿ.'],
            [6.5, 'Cut the torus along imaginary time, between two sites: one leg is cut. What is left is the matrix Xⁿ with two open legs. Divide by Zₙ so its trace is 1: σ = Xⁿ / Zₙ, with eigenvalues pᵢ = λᵢⁿ / Zₙ.'],
            [15.5, 'The purity of σ joins two copies leg to leg and closes the ring: Tr σ² = Tr X²ⁿ / Zₙ² = Z₂ₙ / Zₙ², a ring twice as long. This is the spatial purity Pₛ(n, β).'],
            [23.5, 'Cut along space instead, at one moment τ, through every site. What is left is e^(−βH), with open legs above and below. Divide by Zₙ(β): the thermal state ρ, whose eigenvalues are the Boltzmann probabilities.'],
            [32, 'The purity of ρ stacks two copies and closes time: Tr ρ² = Zₙ(2β) / Zₙ(β)², a torus twice as tall. This is the thermal purity Pτ(n, β).'],
            [39, 'Both purities are ratios of partition functions of small tori: numbers a computer can evaluate.']],
    draw(K, t) {
      const {th} = K, c = 22;
      // left: the cut along imaginary time
      const tx = 90, ty = 66, cut1 = seg(t, 7, 7.8);
      K.torus(tx, ty, 6, 3, c, {fill: th.spaceF});
      K.line(tx + 3 * c, ty - 8, tx + 3 * c, ty + 3 * c + 8, {color: th.space, w: 4, alpha: cut1});
      K.text('cut along τ', tx + 3 * c + 10, ty - 12, {color: th.space, size: 14, alpha: cut1});
      const ry = 220, rx = i => 110 + i * 64, r0 = seg(t, .5, 1.3);
      K.path([[400, ry], [415, ry], [415, ry + 32], [65, ry + 32], [65, ry], [80, ry]], {color: th.ink, w: 2, alpha: r0 * (1 - seg(t, 7.5, 8.5))});
      K.line(80, ry, 400, ry, {color: th.ink, w: 2, alpha: r0});
      for (let i = 0; i < 5; i++) K.box(rx(i), ry, 34, 26, 'X', {alpha: r0});
      K.dot(80, ry, 5, th.space, cut1); K.dot(400, ry, 5, th.space, cut1);
      K.tex('Z_n=\\operatorname{Tr}X^n', 240, 290, {align: 'center', size: 22, alpha: seg(t, 2, 2.8) * (1 - seg(t, 7, 7.5))});
      K.tex('\\sigma=\\frac{X^n}{Z_n}', 240, 295, {align: 'center', size: 22, color: th.space, alpha: seg(t, 9.5, 10.3)});
      K.text('Tr σ = 1,   eigenvalues pᵢ = λᵢⁿ / Zₙ', 70, 352, {size: 15, alpha: seg(t, 11, 11.7)});
      const s2 = seg(t, 16, 16.8), qy = 405;
      K.path([[406, qy], [418, qy], [418, qy + 26], [58, qy + 26], [58, qy], [70, qy]], {color: th.ink, w: 2, alpha: s2});
      K.line(70, qy, 406, qy, {color: th.ink, w: 2, alpha: s2});
      for (let i = 0; i < 10; i++) K.box(87 + i * 34, qy, 24, 22, 'X', {size: 12, color: i < 5 ? th.proved : th.space, fill: i < 5 ? th.provedF : th.spaceF, alpha: s2});
      K.text('copy 1', 155, qy - 18, {align: 'center', color: th.proved, size: 13, alpha: s2}); K.text('copy 2', 325, qy - 18, {align: 'center', color: th.space, size: 13, alpha: s2});
      K.tex('\\operatorname{Tr}\\sigma^2=\\frac{Z_{2n}}{Z_n^{\\,2}}=P_{\\mathrm s}(n,\\beta)', 240, 470, {align: 'center', size: 20, color: th.space, alpha: seg(t, 18, 18.8)});
      // right: the cut at one moment of imaginary time
      const ux = 590, uy = 66, cut2 = seg(t, 24, 24.8);
      K.torus(ux, uy, 6, 3, c, {fill: th.timeF});
      K.line(ux - 8, uy + c, ux + 6 * c + 8, uy + c, {color: th.time, w: 4, alpha: cut2});
      K.text('cut at one moment τ', ux + 6 * c + 14, uy + c + 5, {color: th.time, size: 14, alpha: cut2});
      const ex = 720, ey = 220, e0 = seg(t, 24.5, 25.3);
      for (let i = 0; i < 5; i++) { const x = 615 + i * 52; K.line(x, ey - 38, x, ey + 38, {color: th.ink, w: 2, alpha: e0}); K.dot(x, ey - 38, 4, th.time, e0); K.dot(x, ey + 38, 4, th.time, e0); }
      K.box(ex, ey, 250, 40, 'e^(−βH)', {color: th.time, fill: th.timeF, size: 17, alpha: e0});
      K.tex('\\rho=\\frac{e^{-\\beta H}}{Z_n(\\beta)}', ex, 295, {align: 'center', size: 22, color: th.time, alpha: seg(t, 26, 26.8)});
      K.text('Tr ρ = 1,   eigenvalues e^(−βE) / Z', 560, 352, {size: 15, alpha: seg(t, 27.5, 28.2)});
      const s3 = seg(t, 32.5, 33.3);
      for (let i = 0; i < 5; i++) { const x = 625 + i * 48; K.line(x, 376, x, 446, {color: th.ink, w: 2, alpha: s3}); }
      K.box(ex, 397, 230, 22, 'e^(−βH)', {color: th.time, fill: th.timeF, size: 13, alpha: s3});
      K.box(ex, 425, 230, 22, 'e^(−βH)', {color: th.time, fill: th.timeF, size: 13, alpha: s3});
      K.path([[850, 376], [866, 376], [866, 446], [850, 446]], {color: th.time, w: 2, dash: [4, 3], alpha: s3});
      K.text('Tr', 872, 416, {color: th.time, size: 14, alpha: s3});
      K.tex('\\operatorname{Tr}\\rho^2=\\frac{Z_n(2\\beta)}{Z_n(\\beta)^2}=P_\\tau(n,\\beta)', ex, 482, {align: 'center', size: 20, color: th.time, alpha: seg(t, 34, 34.8)});
      K.text('both purities are ratios of partition functions: computable', W / 2, 40, {align: 'center', italic: true, size: 17, alpha: seg(t, 39.5, 40.3)});
    }
  });

  // 7 ------------------------------------------------------------------
  const fsq = z => z * z / (2 * (1 - z) * (1 - z));
  chapter({
    id: 'square', title: 'Squaring sharpens', dur: 27,
    beats: [[0, 'Square every probability and renormalize. The leader pulls ahead.'],
            [5.5, 'Why this is gluing: σ² / Tr σ² = X²ⁿ / Z₂ₙ is the state of the ring twice as long, and its probabilities are pᵢ² / Σ p². Gluing two copies squares every probability (and the same holds in time).'],
            [13, 'A short lemma bounds what squaring does to the defect ε = 1 − P: the new defect is at most f(ε) = ε² / (2(1 − ε)²), about ε²/2.'],
            [19.5, 'Start below ½ and repeated squaring drives the defect to zero: the cobweb drops off the diagonal.']],
    draw(K, t) {
      const {th} = K;
      let p = [.5, .3, .2]; const k = seg(t, .8, 18) * 3, n = Math.floor(k), fr = ease(k - n);
      const sq = q => { const s = q.reduce((u, v) => u + v * v, 0); return q.map(v => v * v / s); };
      for (let i = 0; i < Math.min(n, 3); i++) p = sq(p);
      const nx = n < 3 ? sq(p) : p, cur = p.map((v, i) => lerp(v, nx[i], n < 3 ? fr : 0));
      K.bars(70, 110, 300, 300, cur, th.ink, {labels: true});
      const P = cur.reduce((u, v) => u + v * v, 0);
      K.text(`squarings: ${Math.min(n, 3)}`, 70, 450, {size: 16});
      K.text(`defect 1 − P = ${(1 - P).toFixed(4)}`, 70, 478, {size: 16, color: th.proved});
      const wy = 1 - seg(t, 12.4, 13);
      K.tex('\\frac{\\sigma^2}{\\operatorname{Tr}\\sigma^2}=\\frac{X^{2n}}{Z_{2n}}', 690, 170, {align: 'center', size: 26, alpha: wy * seg(t, 6, 6.8)});
      K.text('the state of the ring twice as long', 520, 240, {size: 17, alpha: wy * seg(t, 7, 7.7)});
      K.text('its probabilities:  pᵢ  →  pᵢ² / Σ p²', 520, 274, {size: 17, alpha: wy * seg(t, 8.5, 9.2)});
      K.text('two copies glued = one squaring', 520, 316, {size: 18, color: th.proved, alpha: wy * seg(t, 10, 10.7)});
      const gx = 480, gy = 100, gw = 410, gh = 330, zm = .5, X = z => gx + z / zm * gw, Y = z => gy + gh - z / zm * gh, b = seg(t, 13, 13.8);
      K.axes(gx, gy, gw, gh, {alpha: b});
      K.line(X(0), Y(0), X(zm), Y(zm), {color: th.mute, dash: [4, 4], alpha: b, w: 1.5});
      const fp = []; for (let z = 0; z <= zm + 1e-9; z += .005) fp.push([X(z), Y(fsq(z))]); K.path(fp, {color: th.proved, w: 3, alpha: b});
      K.tex('f(z)=\\frac{z^2}{2(1-z)^2}', gx + 30, gy + 30, {size: 20, color: th.proved, alpha: b});
      K.text('old defect →', gx + gw, gy + gh + 24, {align: 'right', color: th.mute, size: 13, alpha: b});
      K.text('new defect', gx - 8, gy - 8, {color: th.mute, size: 13, alpha: b});
      let z = .45; const pts = [[X(z), Y(0)]]; for (let i = 0; i < 5; i++) { const fz = fsq(z); pts.push([X(z), Y(fz)], [X(fz), Y(fz)]); z = fz; }
      const cw = seg(t, 19.5, 25) * (pts.length - 1), ni = Math.floor(cw), sh = pts.slice(0, ni + 1);
      if (ni < pts.length - 1) sh.push([lerp(pts[ni][0], pts[ni + 1][0], cw - ni), lerp(pts[ni][1], pts[ni + 1][1], cw - ni)]);
      K.path(sh, {color: th.ink, w: 2, alpha: b});
    }
  });

  // 8 ------------------------------------------------------------------
  chapter({
    id: 'couple', title: 'Space and time are tied', dur: 31,
    beats: [[0, 'Squaring grows the torus in one direction at a time. We need both.'],
            [4, 'Bookkeeping: every gluing multiplies by a purity. Side by side, Z₂ₙ(β) = Pₛ(n, β) · Zₙ(β)². Stacked, Zₙ(2β) = Pτ(n, β) · Zₙ(β)². These are just the two definitions, rearranged.'],
            [11.5, 'Start from four copies, Zₙ(β)⁴, and build one torus twice as long and twice as tall. Rows first: glue side by side (× Pₛ(n, β)²), then stack (× Pτ(2n, β)).'],
            [18.5, '… or columns first: stack (× Pτ(n, β)²), then glue side by side (× Pₛ(n, 2β)). Both routes end at the same torus, Z₂ₙ(2β).'],
            [25, 'So the two products are equal: Pₛ(n, 2β) Pτ(n, β)² = Pₛ(n, β)² Pτ(2n, β). This coupling identity ties the spatial and thermal purities together.']],
    draw(K, t) {
      const {th} = K, c = 14, L = 150, R = 600, T = 70, B = 290;
      const g1 = ease(seg(t, 12.5, 14)), g2 = ease(seg(t, 14.5, 16.5)), h1 = ease(seg(t, 19, 20.5)), h2 = ease(seg(t, 21, 22.5));
      const tiles = (x, y, cols, rows, nx, ny, al, fill) => { for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) K.torus(x + i * (cols * c + 14), y + j * (rows * c + 14), cols, rows, c, {alpha: al, fill}); };
      tiles(L, T, 4, 3, 2, 2, 1);
      tiles(R, T, 8, 3, 1, 2, Math.max(.2, g1));
      tiles(L, B, 4, 6, 2, 1, Math.max(.2, h1));
      K.torus(R, B, 8, 6, c, {alpha: Math.max(.2, g2, h2), fill: th.provedF});
      K.text('Zₙ(β)⁴', L, T - 14, {color: th.mute, size: 16, alpha: seg(t, 12, 12.7)});
      K.text('Z₂ₙ(2β)', R, B + 6 * c + 26, {color: th.proved, size: 16, alpha: seg(t, 16.5, 17.2)});
      K.arrow(L + 140, T + 30, R - 20, T + 30, {color: th.space, w: 3, alpha: lerp(.25, 1, g1)});
      K.arrow(R + 56, T + 104, R + 56, B - 12, {color: th.time, w: 3, alpha: lerp(.25, 1, g2)});
      K.arrow(L + 56, T + 104, L + 56, B - 12, {color: th.time, w: 3, alpha: lerp(.25, 1, h1)});
      K.arrow(L + 140, B + 42, R - 20, B + 42, {color: th.space, w: 3, alpha: lerp(.25, 1, h2)});
      K.text('× Pₛ(n,β)²', (L + R) / 2 + 70, T + 20, {align: 'center', color: th.space, size: 15, alpha: lerp(.25, 1, g1)});
      K.text('× Pτ(2n,β)', R + 70, (T + B) / 2 + 50, {color: th.time, size: 15, alpha: lerp(.25, 1, g2)});
      K.text('× Pτ(n,β)²', L + 42, (T + B) / 2 + 50, {align: 'right', color: th.time, size: 15, alpha: lerp(.25, 1, h1)});
      K.text('× Pₛ(n,2β)', (L + R) / 2 + 70, B + 32, {align: 'center', color: th.space, size: 15, alpha: lerp(.25, 1, h2)});
      const bk = 1 - seg(t, 11, 11.5);
      K.tex('Z_{2n}(\\beta)=P_{\\mathrm s}(n,\\beta)\\,Z_n(\\beta)^2', W / 2, 440, {align: 'center', size: 21, color: th.space, alpha: bk * seg(t, 5, 5.8)});
      K.tex('Z_{n}(2\\beta)=P_\\tau(n,\\beta)\\,Z_n(\\beta)^2', W / 2, 484, {align: 'center', size: 21, color: th.time, alpha: bk * seg(t, 7, 7.8)});
      K.tex('P_{\\mathrm s}(n,2\\beta)\\,P_\\tau(n,\\beta)^2=P_{\\mathrm s}(n,\\beta)^2\\,P_\\tau(2n,\\beta)', W / 2, 470, {align: 'center', size: 22, alpha: seg(t, 25.5, 26.5)});
    }
  });

  // 9 ------------------------------------------------------------------
  const RPM = 700, RC = [560, 610, 660], REV = [[0, 1, 140], [0, 1, 250], [1, 2, 170], [1, 2, 300], [0, 1, 330]], RX = [190, 230, 280];
  function rpHalf(K, s, col, colF, alpha, lab) {
    const X = x => RPM + (x - RPM) * s;
    K.rect(Math.min(X(RC[0] - 24), X(RPM - 3)), 110, Math.abs(X(RPM - 3) - X(RC[0] - 24)), 250, {fill: colF, alpha});
    for (const c of RC) K.line(X(c), 110, X(c), 360, {color: K.th.mute, w: 1.5, alpha});
    for (const [i, j, y] of REV) { K.line(X(RC[i]), y, X(RC[j]), y, {color: col, w: 3, alpha}); K.dot(X(RC[i]), y, 4, col, alpha); K.dot(X(RC[j]), y, 4, col, alpha); }
    for (const y of RX) { K.line(X(RC[2]), y, X(RPM), y, {color: col, w: 3, alpha}); K.dot(X(RC[2]), y, 4, col, alpha); if (lab) K.text(lab, X(RC[2]) + (s > 0 ? -10 : 10), y + 5, {color: col, size: 14, align: s > 0 ? 'right' : 'left', alpha}); }
  }
  // reflection positivity, rebuilt around worked numbers
  const SYMX = [[1, 1.4, .6327], [2, 1.24, .8314], [3, 1.232, .9385], [4, 1.2784, .9794], [5, 1.3462, .9933], [6, 1.4245, .9979], [7, 1.5097, .9993], [8, 1.6007, .9998]];
  const ROTX = [[1, 1.8, .4012], [2, 1.3, -.1135], [3, .594, -4.1734], [4, -.1918, -47.13], [5, -.9214, -1.0229], [6, -1.4725, .2317], [7, -1.7567, .4764], [8, -1.7338, .4786]];
  chapter({
    id: 'rpwhy', title: 'Why the proof needs a symmetric transfer matrix', dur: 34,
    beats: [[0, 'The squaring step reads pᵢ = λᵢⁿ / Zₙ as probabilities. Let us test that on a 2×2 transfer matrix, with actual numbers.'],
            [6, 'A symmetric matrix. Its eigenvalues, 1.06 and 0.34, are real. Compute Zₙ = Tr Xⁿ and the spatial purity Z₂ₙ/Zₙ² for n = 1 … 8.'],
            [14, 'Every purity lies between 0 and 1 and climbs to 1. This is the behavior the bootstrap relies on.'],
            [19, 'Now a matrix that is not symmetric: a small rotation. Its eigenvalues 0.9 ± 0.4i are complex.'],
            [24.5, 'Zₙ changes sign, and the "purity" becomes negative or huge (−47 at n = 4). The pᵢ are no longer probabilities, so squaring proves nothing. The proof needs X symmetric.']],
    draw(K, t) {
      const {th} = K;
      const table = (x, rows, at, col) => {
        K.text('n', x, 150, {color: th.mute, size: 14, italic: true}); K.text('Zₙ', x + 50, 150, {color: th.mute, size: 14, italic: true}); K.text('Z₂ₙ / Zₙ²', x + 140, 150, {color: th.mute, size: 14, italic: true});
        K.line(x, 158, x + 230, 158, {color: th.line, w: 1});
        rows.forEach((r, i) => {
          const al = seg(t, at + i * .55, at + .4 + i * .55), y = 182 + i * 26, bad = r[2] < 0 || r[2] > 1;
          K.text(String(r[0]), x, y, {size: 16, alpha: al}); K.text(r[1].toFixed(4), x + 50, y, {size: 16, alpha: al, color: r[1] < 0 ? th.fail : th.ink});
          K.text(Math.abs(r[2]) > 10 ? r[2].toFixed(1) : r[2].toFixed(4), x + 140, y, {size: 16, alpha: al, color: bad ? th.fail : col});
        });
      };
      // left: symmetric
      const a = seg(t, 6, 6.8);
      K.text('symmetric', 120, 70, {size: 20, color: th.ok, alpha: a});
      K.text('X = [[0.9, 0.3], [0.3, 0.5]]', 120, 100, {size: 16, alpha: a});
      K.text('eigenvalues 1.06, 0.34 (real)', 120, 124, {size: 14, color: th.mute, alpha: a});
      if (t >= 7) table(120, SYMX, 7, th.ok);
      K.text('all between 0 and 1  ✓', 120, 410, {size: 17, color: th.ok, alpha: seg(t, 14, 15)});
      // right: rotation
      const b = seg(t, 19, 19.8);
      K.text('not symmetric', 560, 70, {size: 20, color: th.fail, alpha: b});
      K.text('X = [[0.9, −0.4], [0.4, 0.9]]', 560, 100, {size: 16, alpha: b});
      K.text('eigenvalues 0.9 ± 0.4i (complex)', 560, 124, {size: 14, color: th.mute, alpha: b});
      if (t >= 20) table(560, ROTX, 20, th.ink);
      K.text('negative, or far above 1  ✗', 560, 410, {size: 17, color: th.fail, alpha: seg(t, 24.5, 25.5)});
      K.text('the purities of a genuine state are always between 0 and 1', W / 2, 470, {align: 'center', size: 17, italic: true, alpha: seg(t, 28, 29)});
      if (t < 6) {
        K.text('pᵢ = λᵢⁿ / Zₙ', W / 2, 230, {align: 'center', size: 34, alpha: seg(t, .5, 1.5)});
        K.text('must be probabilities: real, ≥ 0, adding to 1', W / 2, 280, {align: 'center', size: 18, italic: true, alpha: seg(t, 2, 3)});
      }
    }
  });

  chapter({
    id: 'rpwhere', title: 'Where the symmetry comes from', dur: 47,
    beats: [[0, 'In the proof X is a matrix: its row index is its left leg (the history a of the bond on the site’s left) and its column index is its right leg (the history b of the bond on its right). Symmetric means X(a, b) = X(b, a): mirroring the box left ↔ right changes nothing.'],
            [8, 'Test the mirror on one site. An event’s side is the bond it arrives through: the left bond, shared with the left neighbor, or the right bond, shared with the right neighbor. Example: a = {G_z at 0.2, G_x at 0.6} on the left, b = {G_y at 0.4} on the right.'],
            [15.5, 'For Heisenberg an event puts the same matrix G_α on both sites of its bond, so the site gets G_α whichever side the event comes from. In time order it sees G_z, G_y, G_x: X(a, b) = Tr(G_x G_y G_z) = 1.'],
            [22, 'Mirror it: the left events move to the right and the right event to the left. The site sees the same three matrices at the same times, so X(b, a) = 1 too. Heisenberg’s X is symmetric.'],
            [27.5, 'Counter-example: a bond with −h = Q ⊗ P, where P = diag(1, 0, 0) projects onto |x⟩ and Q = diag(1, 1, 0) = −G_z². A site is the right end of its left bond, so it receives P from there, and the left end of its right bond, so it receives Q from there.'],
            [36.5, 'Two events through the left bond give Tr(P P) = 1; the same two events through the right bond give Tr(Q Q) = 2. So X(a, ∅) ≠ X(∅, a): not symmetric. Reflection positivity is the condition that rules this out: the bond must look the same from both sides, with positive weights.']],
    draw(K, t) {
      const {th} = K, sx = 300, y0 = 110, y1 = 400, Yt = f => y1 - f * (y1 - y0);
      // intro: X as a box whose two legs are its row and column indices
      const I = 1 - seg(t, 7.4, 8);
      if (I > 0) {
        const xbox = (x, l, r, lc, rc, lab, al) => {
          K.line(x - 80, 230, x + 80, 230, {color: th.ink, w: 2, alpha: al});
          K.line(x - 80, 230, x, 230, {color: lc, w: 3, alpha: al}); K.line(x, 230, x + 80, 230, {color: rc, w: 3, alpha: al});
          K.box(x, 230, 50, 40, 'X', {alpha: al});
          K.text(l, x - 88, 236, {align: 'right', color: lc, size: 20, alpha: al}); K.text(r, x + 88, 236, {color: rc, size: 20, alpha: al});
          K.text(lab, x, 290, {align: 'center', size: 20, alpha: al});
        };
        K.text('row index = left leg,   column index = right leg', W / 2, 130, {align: 'center', color: th.mute, size: 17, alpha: I * seg(t, 1, 1.7)});
        xbox(280, 'a', 'b', th.proved, th.time, 'X(a, b)', I * seg(t, .3, 1));
        K.text('mirror', W / 2, 222, {align: 'center', color: th.mute, size: 16, alpha: I * seg(t, 3.5, 4.2)});
        K.text('↔', W / 2, 250, {align: 'center', color: th.mute, size: 26, alpha: I * seg(t, 3.5, 4.2)});
        xbox(680, 'b', 'a', th.time, th.proved, 'X(b, a)', I * seg(t, 3.5, 4.2));
        K.text('symmetric:   X(a, b) = X(b, a)   for all histories a, b', W / 2, 370, {align: 'center', color: th.ok, size: 19, alpha: I * seg(t, 5, 5.7)});
      }
      // one site, its two neighbors, and the events that reach it
      const A = seg(t, 7.6, 8.2) * (1 - seg(t, 27, 27.5));
      if (A > 0) {
        const sw = ease(seg(t, 22.5, 24.5));
        [-1, 1].forEach(s => { K.line(sx + s * 140, y0, sx + s * 140, y1, {color: th.mute, w: 1.5, alpha: A * .6});
          K.text(s < 0 ? 'left neighbor' : 'right neighbor', sx + s * 140, y1 + 26, {align: 'center', color: th.mute, size: 13, alpha: A}); });
        K.line(sx, y0, sx, y1, {color: th.ink, w: 3, alpha: A});
        K.text('the site', sx, y1 + 26, {align: 'center', size: 14, alpha: A});
        K.text('left bond', sx - 70, y0 - 16, {align: 'center', color: th.proved, size: 15, alpha: A});
        K.text('right bond', sx + 70, y0 - 16, {align: 'center', color: th.time, size: 15, alpha: A});
        K.text('τ ↑', sx - 260, y0 + 6, {color: th.time, size: 14, alpha: A});
        [[.2, -1, 'z'], [.6, -1, 'x'], [.4, 1, 'y']].forEach(([f, s0, g], i) => {
          const s = lerp(s0, -s0, sw), y = Yt(f), al = A * seg(t, 8.5 + i * .8, 9.2 + i * .8), col = s0 < 0 ? th.proved : th.time, both = A * seg(t, 16, 16.7);
          K.line(sx + s * 140, y, sx, y, {color: col, w: 3, alpha: al}); K.dot(sx + s * 140, y, 5, col, al); K.dot(sx, y, 5, th.ink, al);
          K.text(`G_${g}`, sx + (s > 0 ? 10 : -10), y - 8, {align: s > 0 ? 'left' : 'right', size: 16, alpha: al});
          K.text(`G_${g}`, sx + s * 140 + (s > 0 ? -10 : 10), y - 8, {align: s > 0 ? 'right' : 'left', size: 16, color: th.mute, alpha: Math.min(al, both)});
          K.text(`τ = ${f}`, sx + s * 150 + (s > 0 ? 6 : -6), y + 5, {align: s > 0 ? 'left' : 'right', color: th.mute, size: 13, alpha: al});
        });
        K.text('same G at both ends of an event', 560, 130, {size: 16, color: th.mute, alpha: A * seg(t, 16, 16.7)});
        K.text('time-ordered product (later on the left):', 560, 170, {size: 16, alpha: A * seg(t, 17, 17.7)});
        K.text('X(a, b) = Tr(G_x G_y G_z) = 1', 560, 205, {size: 22, alpha: A * seg(t, 18, 18.8)});
        K.text('sides swapped: same matrices, same times', 560, 260, {size: 16, color: th.ok, alpha: A * seg(t, 24.5, 25.2)});
        K.text('X(b, a) = 1 = X(a, b)', 560, 295, {size: 22, color: th.ok, alpha: A * seg(t, 25, 25.8)});
      }
      // a bond that is not mirror symmetric
      const B = seg(t, 27.5, 28.3);
      if (B > 0) {
        K.text('a bond that is not the same from both sides:', 80, 80, {size: 17, alpha: B});
        const bb = B * seg(t, 28.5, 29.2);
        K.line(130, 130, 250, 130, {color: th.fail, w: 3, alpha: bb}); K.dot(130, 130, 6, th.ink, bb); K.dot(250, 130, 6, th.ink, bb);
        K.text('Q', 130, 114, {align: 'center', size: 16, alpha: bb}); K.text('P', 250, 114, {align: 'center', size: 16, alpha: bb});
        K.text('left site', 130, 154, {align: 'center', color: th.mute, size: 13, alpha: bb}); K.text('right site', 250, 154, {align: 'center', color: th.mute, size: 13, alpha: bb});
        K.text('its event: Q on the left end, P on the right end', 300, 126, {size: 16, color: th.fail, alpha: bb});
        K.text('P = diag(1, 0, 0) projects onto |x⟩;   Q = diag(1, 1, 0) = −G_z²', 300, 152, {size: 15, color: th.mute, alpha: B * seg(t, 30, 30.7)});
        const panel = (x0, side, nb, me, val, at) => {
          const al = B * seg(t, at, at + .8), x = x0 + 120;
          K.line(x + side * 120, 210, x + side * 120, 400, {color: th.mute, w: 1.5, alpha: al * .6});
          K.line(x, 210, x, 400, {color: th.ink, w: 3, alpha: al});
          [.3, .7].forEach(f => { const y = 400 - f * 190; K.line(x + side * 120, y, x, y, {color: side < 0 ? th.proved : th.time, w: 3, alpha: al}); K.dot(x, y, 5, th.ink, al); K.dot(x + side * 120, y, 4, th.mute, al);
            K.text(me, x + (side > 0 ? 10 : -10), y - 8, {align: side > 0 ? 'left' : 'right', size: 16, alpha: al});
            K.text(nb, x + side * 120 + (side > 0 ? -10 : 10), y - 8, {align: side > 0 ? 'right' : 'left', color: th.mute, size: 15, alpha: al}); });
          K.text(side < 0 ? 'two events through the left bond' : 'the same events through the right bond', x, 428, {align: 'center', color: th.mute, size: 14, alpha: al});
          K.text(val, x, 466, {align: 'center', size: 22, alpha: B * seg(t, 37, 37.8) * (side < 0 ? 1 : seg(t, 38, 38.8))});
        };
        panel(140, -1, 'Q', 'P', 'Tr(P P) = 1', 31.5);
        panel(520, 1, 'P', 'Q', 'Tr(Q Q) = 2', 33.5);
        K.text('X(a, ∅) = 1  ≠  X(∅, a) = 2', W / 2, 508, {align: 'center', size: 18, color: th.fail, alpha: seg(t, 39.5, 40.3)});
      }
    }
  });

  // realignment test for a spin-1 bond
  const G3 = [0, 1, 2].map(a => { const g = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; const b = (a + 1) % 3, c = (a + 2) % 3; g[b][c] = 1; g[c][b] = -1; return g; });
  function minusH(sign, dz = 1) { // -h = sign * sum_a c_a G_a (x) G_a; sign +1 is the antiferromagnet, since S (x) S = -G (x) G
    const M = Array.from({length: 9}, () => new Array(9).fill(0));
    for (let a = 0; a < 3; a++) { const w = sign * (a === 2 ? dz : 1);
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) M[i * 3 + j][k * 3 + l] += w * G3[a][i][k] * G3[a][j][l]; }
    return M;
  }
  function realign(M) { const R = Array.from({length: 9}, () => new Array(9).fill(0)); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) R[i * 3 + k][j * 3 + l] = M[i * 3 + j][k * 3 + l]; return R; }
  function eigSym(A0) {
    const n = A0.length, A = A0.map(r => r.slice());
    for (let sw = 0; sw < 60; sw++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j]; if (off < 1e-20) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-15) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]), tt = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(tt * tt + 1), s = tt * c;
        for (let k = 0; k < n; k++) { const a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
      }
    }
    return A.map((r, i) => r[i]).sort((x, y) => y - x);
  }
  const IX = ['x', 'y', 'z'];
  // Four legs of a two-site tensor: i, k (out, in) on the left site, j, l on the right. m = 0 points
  // them up (outputs) and down (inputs); m = 1 bends them to the left and right, the realigned reading.
  function bondLegs(K, cx, cy, hw, m, al) {
    const {th} = K;
    [['i', -1, -1], ['k', -1, 1], ['j', 1, -1], ['l', 1, 1]].forEach(([nm, sd, ud]) => {
      const ax = cx + sd * hw, ay = cy + ud * 15, ex = lerp(ax, cx + sd * (hw + 60), m), ey = lerp(cy + ud * 50, cy + ud * 15, m);
      const dx = ex - ax, dy = ey - ay, d = Math.hypot(dx, dy) || 1, col = sd < 0 ? th.proved : th.time;
      K.line(ax, ay, ex, ey, {color: col, w: 2.5, alpha: al});
      K.text(nm, ex + dx / d * 12, ey + dy / d * 12, {align: 'center', base: 'middle', color: col, size: 16, italic: true, alpha: al});
    });
  }
  chapter({
    id: 'rptest', title: 'Testing a bond for reflection positivity', dur: 45,
    beats: [[0, 'A test any bond can be put through. For Heisenberg the bond is h = S_j·S_{j+1}, so −h = Σ_α G_α ⊗ G_α. As a tensor, −h has four legs: i, k on the left site (out, in) and j, l on the right site. Its entries are ⟨i j| −h |k l⟩ = Σ_α (G_α)_{ik} (G_α)_{jl}.'],
            [8.5, 'The usual matrix reading: rows (i, j) are the outputs, columns (k, l) the inputs. That is this 9 × 9 matrix, −h for two spin-1 sites. Example: ⟨x x| −h |y y⟩ = (G_z)_{xy} (G_z)_{xy} = 1.'],
            [16, 'Realignment: read the same 81 numbers with the legs grouped differently, rows (i, k) from the left site and columns (j, l) from the right. Nothing is computed; each entry just moves: ⟨x x| −h |y y⟩ goes to row (x, y), column (x, y).'],
            [24.5, 'Why this grouping: −h is two boxes G_α joined by a summed leg α. Split between the sites, each half is a vector v_α = vec(G_α), and the joined leg is a matrix product: realigned −h = Σ_α v_α v_αᵀ, a sum of squares.'],
            [33, 'A sum of squares has no negative eigenvalues. Here they are 2, 2, 2, 0, …: the antiferromagnet passes.'],
            [38, 'Ferromagnet: h = −S·S, so −h = −Σ G_α ⊗ G_α and the realigned matrix has eigenvalues −2: not a sum of squares, so it fails the test in this basis. In general: realigned −h ⪰ 0, up to an on-site shift, means reflection positive.']],
    draw(K, t) {
      const {th} = K, fm = t >= 38, M = minusH(fm ? -1 : 1), cs = 30, ox = 140, oy = 80;
      const mv = ease(seg(t, 20.5, 23)), one = ease(seg(t, 19, 20.5)), hl = t >= 10 && t < 24.5;
      K.text(fm ? 'ferromagnet:  h = −S_j·S_{j+1}' : 'Heisenberg antiferromagnet:  h = S_j·S_{j+1}', 480, 62, {color: fm ? th.fail : th.mute, size: 16, alpha: seg(t, .3, 1)});
      // the 9 x 9 matrix, before and after realignment
      K.ctx.save(); K.ctx.globalAlpha *= seg(t, 8.5, 9.2);
      for (let r = 0; r < 9; r++) {
        const p = Math.floor(r / 3), q = r % 3;
        const rowL = mv < .5 ? [th.proved, th.time] : [th.proved, th.proved], colL = mv < .5 ? [th.proved, th.time] : [th.time, th.time];
        K.text(IX[p], ox - 26, oy + r * cs + cs / 2 + 5, {align: 'center', color: rowL[0], size: 12});
        K.text(IX[q], ox - 14, oy + r * cs + cs / 2 + 5, {align: 'center', color: rowL[1], size: 12});
        K.text(IX[p], ox + r * cs + cs / 2 - 5, oy - 8, {align: 'center', color: colL[0], size: 12});
        K.text(IX[q], ox + r * cs + cs / 2 + 5, oy - 8, {align: 'center', color: colL[1], size: 12});
      }
      K.text(mv < .5 ? 'rows (i, j)   columns (k, l)' : 'rows (i, k)   columns (j, l)', ox, oy + 9 * cs + 28, {color: th.mute, size: 14});
      K.rect(ox - 2, oy - 2, 9 * cs + 4, 9 * cs + 4, {stroke: th.line, w: 1.5});
      for (let a = 1; a < 3; a++) { K.line(ox + a * 3 * cs, oy, ox + a * 3 * cs, oy + 9 * cs, {color: th.line, w: 1}); K.line(ox, oy + a * 3 * cs, ox + 9 * cs, oy + a * 3 * cs, {color: th.line, w: 1}); }
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) {
        const v = M[i * 3 + j][k * 3 + l]; if (Math.abs(v) < 1e-9) continue;
        const special = i === 0 && j === 0 && k === 1 && l === 1, f = special ? Math.max(one, mv) : mv;
        const r = lerp(i * 3 + j, i * 3 + k, f), c = lerp(k * 3 + l, j * 3 + l, f);
        K.rect(ox + c * cs + 2, oy + r * cs + 2, cs - 4, cs - 4, {fill: v > 0 ? th.proved : th.time, alpha: special && hl ? 1 : .8});
        if (special && hl) K.rect(ox + c * cs, oy + r * cs, cs, cs, {stroke: th.ink, w: 2.5});
        K.text(v > 0 ? '+' : '−', ox + c * cs + cs / 2, oy + r * cs + cs / 2 + 5, {align: 'center', size: 13, color: th.bg});
      }
      K.ctx.restore();
      // the tensor -h and its two readings
      const P0 = seg(t, .3, 1) * (1 - seg(t, 24, 24.5)), m = ease(seg(t, 16.5, 18)), cx = 720, cy = 190;
      if (P0 > 0) {
        bondLegs(K, cx, cy, 50, m, P0);
        K.box(cx, cy, 140, 30, '−h', {alpha: P0});
        K.text('left site', cx - 50, cy - 92, {align: 'center', color: th.proved, size: 13, alpha: P0 * seg(t, 1, 1.7) * (1 - seg(t, 8.5, 9))});
        K.text('right site', cx + 50, cy - 92, {align: 'center', color: th.time, size: 13, alpha: P0 * seg(t, 1, 1.7) * (1 - seg(t, 8.5, 9))});
        K.text('rows (i, j)', cx, cy - 82, {align: 'center', size: 15, alpha: P0 * seg(t, 9, 9.7) * (1 - m)});
        K.text('columns (k, l)', cx, cy + 92, {align: 'center', size: 15, alpha: P0 * seg(t, 9, 9.7) * (1 - m)});
        K.text('rows (i, k)', cx - 110, cy + 52, {align: 'center', color: th.proved, size: 15, alpha: P0 * m});
        K.text('columns (j, l)', cx + 110, cy + 52, {align: 'center', color: th.time, size: 15, alpha: P0 * m});
      }
      // -h as two G boxes joined by the summed leg alpha
      const F0 = seg(t, 4.5, 5.2) * (1 - seg(t, 8.2, 8.5)), F3 = seg(t, 24.5, 25.2) * (1 - seg(t, 32.6, 33.2));
      const fact = (y, mm, al) => {
        bondLegs(K, cx, y, 50, mm, al);
        K.line(cx - 50, y, cx + 50, y, {color: th.ink, w: 2.5, alpha: al});
        K.box(cx - 50, y, 46, 30, 'G_α', {size: 14, alpha: al}); K.box(cx + 50, y, 46, 30, 'G_α', {size: 14, alpha: al});
        K.text('α', cx, y - 8, {align: 'center', color: th.ink, size: 15, italic: true, alpha: al});
      };
      if (F0 > 0) {
        fact(360, 0, F0);
        K.text('−h  =  Σ_α', cx - 110, 366, {align: 'right', size: 18, alpha: F0});
        K.tex('\\langle ij|-h|kl\\rangle=\\sum_\\alpha (G_\\alpha)_{ik}(G_\\alpha)_{jl}', cx, 455, {align: 'center', size: 18, alpha: F0 * seg(t, 6, 6.7)});
      }
      const tx = (s, y, a0, a1, o = {}) => K.text(s, 480, y, {size: 16, ...o, alpha: seg(t, a0, a0 + .7) * (1 - seg(t, a1 - .4, a1))});
      tx('⟨x x| −h |y y⟩ = (G_z)_xy · (G_z)_xy = 1', 360, 10, 16);
      tx('every entry: Σ_α (G_α)_ik (G_α)_jl', 392, 11.5, 16, {color: th.mute});
      tx('same tensor, same 81 numbers: only the grouping of legs changes', 360, 17, 24.5);
      tx('⟨x x| −h |y y⟩ = +1   →   row (x, y), column (x, y)', 392, 18.5, 24.5, {color: th.proved});
      K.text('realignment:', 480, 440, {size: 16, alpha: seg(t, 21.5, 22.2) * (1 - seg(t, 24, 24.5))});
      K.tex('R[-h]_{(ik),(jl)}=\\langle ij|{-h}|kl\\rangle', 720, 440, {align: 'center', size: 20, alpha: seg(t, 21.5, 22.2) * (1 - seg(t, 24, 24.5))});
      if (F3 > 0) {
        fact(170, ease(seg(t, 25, 26.3)), F3);
        K.text('v_α = vec(G_α)', cx - 160, 236, {align: 'center', color: th.proved, size: 14, alpha: F3 * seg(t, 26.3, 27)});
        K.text('v_α = vec(G_α)', cx + 160, 236, {align: 'center', color: th.time, size: 14, alpha: F3 * seg(t, 26.3, 27)});
        K.text('joined leg α: a sum over α, i.e. a matrix product', 480, 268, {size: 15, color: th.mute, alpha: F3 * seg(t, 27, 27.7)});
      }
      const D = seg(t, 27.5, 28.3) * (1 - seg(t, 32.6, 33.2));
      if (D > 0) {
        K.text('realigned −h  =', 470, 300, {size: 16, alpha: D});
        for (let a = 0; a < 3; a++) {
          const v = []; for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) v.push(G3[a][i][k]);
          const bx = 470 + a * 150, by = 315, c2 = 12;
          for (let r = 0; r < 9; r++) for (let cc = 0; cc < 9; cc++) { const x = v[r] * v[cc] * (fm ? -1 : 1); if (!x) continue; K.rect(bx + cc * c2, by + r * c2, c2 - 1, c2 - 1, {fill: x > 0 ? th.proved : th.time, alpha: D * seg(t, 28 + a * .6, 28.6 + a * .6)}); }
          K.rect(bx - 1, by - 1, 9 * c2 + 1, 9 * c2 + 1, {stroke: th.line, w: 1, alpha: D});
          K.text(`v_${IX[a]} v_${IX[a]}ᵀ`, bx + 54, by + 9 * c2 + 20, {align: 'center', size: 13, color: th.mute, alpha: D});
          if (a < 2) K.text('+', bx + 9 * c2 + 17, by + 60, {align: 'center', size: 22, alpha: D});
        }
        K.text('each term is a square (v vᵀ): never negative', 470, 478, {color: th.ok, size: 15, alpha: D * seg(t, 30, 30.7)});
      }
      const E2 = seg(t, 33, 33.8);
      if (E2 > 0) {
        const ev = eigSym(realign(M)), gx = 500, gy = 140, gw = 380, gh = 230, Y = v => gy + gh / 2 - v / 2.4 * gh / 2, bw = gw / 9 - 6;
        K.line(gx, Y(0), gx + gw, Y(0), {color: th.mute, w: 1.5, alpha: E2});
        ev.forEach((v, i) => K.rect(gx + i * (bw + 6), Math.min(Y(0), Y(v)), bw, Math.max(2, Math.abs(Y(v) - Y(0))), {fill: v < -1e-9 ? th.fail : th.ok, alpha: E2}));
        K.text('eigenvalues of the realigned matrix', gx, gy - 14, {color: th.mute, size: 14, alpha: E2});
        K.text(fm ? 'ferromagnet: −2, −2, −2  ✗ fails the test in this basis' : 'antiferromagnet: 2, 2, 2, 0, …  ✓ reflection positive', gx, gy + gh + 34, {color: fm ? th.fail : th.ok, size: 16, alpha: E2});
      }
    }
  });

  // 10 -----------------------------------------------------------------
  chapter({
    id: 'double', title: 'One doubling step', dur: 36,
    beats: [[0, 'Start with both purities at least 1 − ε at size (n, β): the spatial purity Pₛ(n, β) and the thermal purity Pτ(2n, β).'],
            [5, 'Move 1, couple. The identity gives Pₛ(n, 2β) = Pₛ(n, β)² Pτ(2n, β) / Pτ(n, β)². Purities are at most 1, so the division only helps: Pₛ(n, 2β) ≥ (1 − ε)³, a defect of about 3ε.'],
            [13, 'Move 2, square in space. Gluing two rings sends the defect 3ε to about (3ε)²/2 = 4.5ε², below 5ε². That bounds Pₛ(2n, 2β).'],
            [19.5, 'Move 3, couple again, now at size 2n: Pτ(4n, β) = Pₛ(2n, 2β) Pτ(2n, β)² / Pₛ(2n, β)² ≥ (1 − 5ε²)(1 − ε)², a defect of about 2ε.'],
            [26.5, 'Move 4, square in time. Stacking sends the defect 2ε to about 2ε², below 5ε². That bounds Pτ(4n, 2β).'],
            [31, 'Out comes the same pair of bounds at (2n, 2β), with ε replaced by 5ε². The step repeats.']],
    draw(K, t) {
      const {th} = K, e = .01, d = .5 * Math.pow(Math.pow(1 - e, -3) - 1, 2);
      const xs = 90, xt = 520, bw = 350, bh = 50, rows = [60, 145, 230, 315, 400];
      const box = (x, r, col, s, at, dash) => { const al = seg(t, at, at + .7); K.rect(x, rows[r], bw, bh, {fill: th.surface, stroke: col, w: 2, dash, alpha: al}); K.text(s, x + 16, rows[r] + 32, {size: 15, alpha: al}); };
      box(xs, 0, th.space, 'Pₛ(n,β) ≥ 1 − ε = 0.99', .3); box(xt, 0, th.time, 'Pτ(2n,β) ≥ 1 − ε = 0.99', .3);
      box(xs, 1, th.space, `Pₛ(n,2β) ≥ (1−ε)³ = ${Math.pow(1 - e, 3).toFixed(4)}`, 9, [6, 5]);
      box(xs, 2, th.space, `Pₛ(2n,2β) ≥ 1 − ${d.toFixed(5)}`, 16);
      box(xt, 3, th.time, `Pτ(4n,β) ≥ ${(Math.pow(1 - e, 2) * (1 - d)).toFixed(4)}`, 23, [6, 5]);
      box(xt, 4, th.time, 'Pτ(4n,2β) ≥ 1 − 5ε² = 0.9995', 28.5);
      const mv = (x, r1, r2, k, lab, at) => { const al = seg(t, at - .8, at); K.arrow(x, rows[r1] + bh, x, rows[r2] - 2, {color: th.proved, w: 2.5, alpha: al}); K.text(`${k}  ${lab}`, x + 12, (rows[r1] + bh + rows[r2]) / 2 + 5, {color: th.proved, size: 14, alpha: al}); };
      mv(xs + 40, 0, 1, '①', 'couple', 9); mv(xs + 40, 1, 2, '②', 'square', 16); mv(xt + 40, 0, 3, '③', 'couple', 23); mv(xt + 40, 3, 4, '④', 'square', 28.5);
      K.arrow(xt, rows[0] + bh / 2, xs + bw - 30, rows[1] + 4, {color: th.mute, w: 1.5, alpha: seg(t, 6, 6.6)});
      K.arrow(xs + bw, rows[2] + bh / 2, xt - 6, rows[3] + 10, {color: th.mute, w: 1.5, alpha: seg(t, 20, 20.6)});
      // the reason for each move, one line at a time
      [[5.5, 12.6, 'Pₛ(n,2β) = Pₛ(n,β)² Pτ(2n,β) / Pτ(n,β)²  ≥  (1 − ε)³      (Pτ ≤ 1)'],
       [13.5, 19, 'squaring lemma:   defect 3ε  →  ≈ (3ε)²/2 = 4.5ε²  ≤  5ε²'],
       [20, 26, 'Pτ(4n,β) = Pₛ(2n,2β) Pτ(2n,β)² / Pₛ(2n,β)²  ≥  (1 − 5ε²)(1 − ε)²'],
       [27, 30.6, 'squaring lemma:   defect 2ε  →  ≈ (2ε)²/2 = 2ε²  ≤  5ε²']].forEach(([a, b, s]) =>
        K.text(s, xs, 495, {size: 17, alpha: seg(t, a, a + .7) * (1 - seg(t, b - .4, b))}));
      K.text('(n, β, ε)  →  (2n, 2β, 5ε²)', xs, 495, {size: 24, color: th.proved, alpha: seg(t, 31.5, 32.5)});
    }
  });

  // 11 -----------------------------------------------------------------
  chapter({
    id: 'iterate', title: 'Repeat, and read off the gap', dur: 25,
    beats: [[0, 'Each step doubles β and roughly squares ε. Plot log ε against β.'],
            [6, 'The points fall on a straight line: ε decays exponentially in β.'],
            [10.5, 'Why that is a gap: a thermal purity ≥ 1 − ε forces the ground state to hold probability p₀ ≥ 1 − ε (purity ≤ largest probability). The first excited level then has p₁ ≤ ε, and p₁/p₀ = e^(−βΔ).'],
            [18, 'So e^(−βΔ) ≤ ε/(1 − ε), and Δ ≥ log((1 − ε)/ε)/β: the slope of the line. For spin 1: Δ(L) > log 20 / 784 ≈ 0.0038 for every even L ≥ 2304.']],
    draw(K, t) {
      const {th} = K, gx = 110, gy = 60, gw = 760, gh = 360, bmax = 13000, lmin = -24;
      const X = b => gx + b / bmax * gw, Y = l => gy + l / lmin * gh;
      K.axes(gx, gy, gw, gh);
      [0, -5, -10, -15, -20].forEach(l => { K.line(gx, Y(l), gx + gw, Y(l), {color: th.line, w: 1}); K.text(l ? `10${String(l).split('').map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(c)]).join('')}` : '1', gx - 8, Y(l) + 5, {align: 'right', color: th.mute, size: 13}); });
      [784, 1568, 3136, 6272, 12544].forEach((b, j) => {
        const l = -Math.pow(2, j) * Math.log10(20) - Math.log10(5), al = seg(t, .5 + j * 1.1, 1.1 + j * 1.1);
        K.dot(X(b), Y(l), 6, th.proved, al);
        K.text(`j=${j}`, X(b), gy + gh + 22, {align: 'center', color: th.mute, size: 13, alpha: al});
      });
      K.text('β →', gx + gw, gy + gh + 44, {align: 'right', color: th.mute, size: 14});
      K.text('defect ε (log scale)', gx + 6, gy - 14, {color: th.mute, size: 14});
      const s = seg(t, 6, 8), slope = -Math.log10(20) / 784;
      K.line(X(0), Y(-Math.log10(5)), X(bmax * s), Y(-Math.log10(5) + slope * bmax * s), {color: th.ink, w: 2, alpha: .8});
      K.text('slope: −log 20 / 784 per unit of β', X(5600), Y(-17.5), {color: th.ink, size: 16, alpha: seg(t, 8, 9)});
      K.text('Pτ ≥ 1 − ε   ⇒   p₀ ≥ 1 − ε', X(4300), Y(-1.4), {size: 17, alpha: seg(t, 11.5, 12.3)});
      K.text('p₁ ≤ ε,    p₁ / p₀ = e^(−βΔ)', X(4300), Y(-3), {size: 17, alpha: seg(t, 14, 14.8)});
      K.text('e^(−βΔ) ≤ ε/(1 − ε)   ⇒   Δ ≥ log((1 − ε)/ε) / β', X(4300), Y(-4.6), {size: 17, color: th.proved, alpha: seg(t, 18.5, 19.3)});
      K.text('Δ(L) > log 20 / 784', X(8200), Y(-8), {color: th.proved, size: 24, alpha: seg(t, 20.5, 21.5)});
    }
  });

  // the seed certificate, run on a toy spectrum with every number shown
  function seedtoyBody(K, t) {
    const {th} = K, lam = [1, .6, -.5, .4];
    const cx = 280, sc = 170, y = 140, X = v => cx + v * sc;
    K.line(X(-1.3), y, X(1.3), y, {color: th.mute, w: 1.5});
    [-1, -.5, 0, .5, 1].forEach(v => K.text(String(v), X(v), y + 30, {align: 'center', color: th.mute, size: 13}));
    lam.forEach((l, i) => { const hid = 1 - seg(t, 41, 42); K.line(X(l), y - 14, X(l), y + 14, {color: i ? th.ink : th.space, w: i ? 2 : 3, alpha: .25 + .75 * seg(t, .5 + i * .4, 1 + i * .4) * (t < 8 ? 1 : .35)}); });
    K.text('hidden eigenvalues', X(-1.3), y - 26, {color: th.mute, size: 13, italic: true});
    const v1 = seg(t, 9, 10);
    if (v1 > 0) { const b = Math.pow(1.066377, 1 / 6); K.rect(X(-b), y - 20, X(b) - X(-b), 40, {fill: th.provedF, stroke: th.proved, w: 1.5, alpha: v1}); K.text('|λ| ≤ 1.011', X(b) + 10, y + 5, {color: th.proved, size: 15, alpha: v1}); }
    const v2 = seg(t, 26.5, 27.5);
    if (v2 > 0) { K.rect(X(.97), y - 22, X(1.05) - X(.97), 44, {fill: th.spaceF, stroke: th.space, w: 1.5, alpha: v2}); K.text('λ₀ > 0.97', X(1.0), y - 34, {align: 'center', color: th.space, size: 15, alpha: v2}); }
    // ledger of numbers
    const rows = [
      [8.5, 'Z₆ = Σ λ⁶', '1.0664', '6-site ring, computed exactly'],
      [10.5, 'max |λ|', '≤ 1.011', 'each λ⁶ ≤ Z₆'],
      [16.5, 'Z₂₄', '≥ 1', '24-site ring: trial-state lower bound'],
      [19, 'if all |λ| ≤ 0.97:  Σ λ²⁴', '≤ 0.484', 'fill the cap 0.97⁶ = 0.833 first'],
      [23, 'contradiction', '⇒ λ₀ > 0.97', ''],
      [27.5, 'rest: Σ λ⁶', '≤ 0.233', '1.0664 − 0.833'],
      [29.5, 'rest: Σ |λ|²⁴', '≤ 0.0030', '0.233⁴'],
      [35, 'Pₛ(24) ≥ (1 + 0.0030/0.97²⁴)⁻²', '> 0.987', 'ε < 0.013'],
    ];
    K.text('ledger', 80, 230, {color: th.mute, size: 14, italic: true});
    K.line(80, 238, 880, 238, {color: th.line, w: 1});
    rows.forEach(([at, a, b, c], i) => { const al = seg(t, at, at + .8), yy = 264 + i * 30;
      K.text(a, 80, yy, {size: 17, alpha: al}); K.text(b, 470, yy, {size: 17, color: i === 4 || i === 7 ? th.space : th.ink, alpha: al}); K.text(c, 600, yy, {size: 14, color: th.mute, italic: true, alpha: al}); });
    // capped-mass picture
    const cm = seg(t, 19.5, 20.5) * (1 - seg(t, 25.6, 26.2));
    if (cm > 0) {
      const bx = 640, by = 70, bh = 120, cap = .833 / 1.0664;
      K.text('sixth-moment mass 1.0664, cap 0.833 each', bx, by - 12, {size: 13, color: th.mute, alpha: cm});
      K.rect(bx, by + bh - cap * bh, 40, cap * bh, {fill: th.proved, alpha: cm * .7}); K.rect(bx + 50, by + bh - (1 - cap) * bh, 40, (1 - cap) * bh, {fill: th.proved, alpha: cm * .7});
      K.line(bx - 6, by + bh - cap * bh, bx + 100, by + bh - cap * bh, {color: th.proved, dash: [4, 3], w: 1.5, alpha: cm});
      K.text('worst case: one full bar, the rest', bx + 104, by + bh - 20, {size: 13, color: th.mute, alpha: cm});
    }
  }
  chapter({
    id: 'seedtoy', title: 'The starting numbers, worked on a toy', dur: 74,
    beats: [[0, 'The doubling step needs a start: a purity close to 1 at one size n₀. The purity is Σ pᵢ² with pᵢ = λᵢⁿ / Zₙ, so it depends on the eigenvalues λ of X.'],
            [7, 'But X acts on whole histories: far too big to diagonalize. What a computer can get exactly is Zₙ = Tr Xⁿ = Σ λⁿ, the partition function of a small ring of n sites. Can a few such numbers pin down the purity?'],
            [14, 'Try it on a toy. Pretend X has four eigenvalues, 1, 0.6, −0.5 and 0.4, and play by the rules: use only the traces Zₙ, never the λ themselves.'],
            [22, 'Step 1, moments. Z₆ = Tr X⁶ = Σλ⁶ is the 6-site ring, small enough to compute exactly: 1.0664. Every term λ⁶ is ≥ 0, so no single |λ| can exceed 1.0664^(1/6) = 1.011.'],
            [30, 'Step 2, a trial state. Z₂₄ = Σλ²⁴ is the 24-site ring, too big to diagonalize. But Z₂₄ = Σₘ e^(−βEₘ) ≥ e^(−βE₀), and any trial state ψ has E₀ ≤ ⟨ψ|H|ψ⟩. If that energy is below zero (energies are shifted per site), then Z₂₄ ≥ 1.'],
            [37, 'Now suppose every |λ| were ≤ 0.97, so each λ⁶ ≤ 0.833, with all λ⁶ adding to 1.0664. Since λ²⁴ = (λ⁶)⁴ and fourth powers reward concentration, the worst case is one λ⁶ at the cap and the rest in one more: Σλ²⁴ ≤ 0.833⁴ + 0.233⁴ = 0.484 < 1. Contradiction.'],
            [47, 'So one eigenvalue λ₀ exceeds 0.97, and λ₀⁶ ≥ 0.833. The others share at most 1.0664 − 0.833 = 0.233 of the sixth moment, and for numbers ≥ 0 a sum of fourth powers is at most the fourth power of the sum: Σ_rest |λ|²⁴ ≤ 0.233⁴ = 0.0030.'],
            [55, 'Step 3, the purity. Recall Pₛ(24) = Z₄₈ / Z₂₄² = Σλ⁴⁸ / (Σλ²⁴)². On top every term is ≥ 0, so keep only λ₀: Σλ⁴⁸ ≥ λ₀⁴⁸. Below, Σλ²⁴ = λ₀²⁴ + rest, and the rest is at most 0.0030.'],
            [62, 'So Pₛ(24) ≥ λ₀⁴⁸ / (λ₀²⁴ + 0.0030)² = (1 + 0.0030/λ₀²⁴)⁻² ≥ (1 + 0.0030/0.97²⁴)⁻² = 0.98765… > 0.987: a certified starting purity, defect ε < 0.013, found without knowing the spectrum (the true value is 0.99999).'],
            [69, 'The real certificate follows the same three steps on the spin-1 chain, with exact integers throughout. The next chapter draws its ingredients.']],
    draw(K, t) {
      const {th} = K, a = 1 - seg(t, 13.4, 14), b = seg(t, 14, 14.6);
      if (a > 0) {
        K.ctx.save(); K.ctx.globalAlpha *= a;
        K.text('what the doubling step needs at the start:', 120, 100, {color: th.mute, size: 16, alpha: seg(t, .3, 1)});
        K.tex('P_{\\mathrm s}=\\sum_i p_i^2,\\qquad p_i=\\frac{\\lambda_i^{\\,n}}{Z_n}\\qquad (\\lambda_i:\\ \\text{eigenvalues of }X)', W / 2, 160, {align: 'center', size: 22, alpha: seg(t, .8, 1.6)});
        K.text('X: rows and columns are whole histories, far too many to diagonalize', W / 2, 240, {align: 'center', color: th.fail, size: 17, alpha: seg(t, 7.3, 8)});
        K.tex('Z_n=\\operatorname{Tr}X^n=\\sum_i\\lambda_i^{\\,n}=\\operatorname{Tr}e^{-\\beta H}\\ \\text{on a ring of }n\\text{ sites}', W / 2, 305, {align: 'center', size: 22, alpha: seg(t, 8.5, 9.3)});
        K.text('computable exactly for small rings', W / 2, 355, {align: 'center', color: th.ok, size: 17, alpha: seg(t, 9.5, 10.2)});
        K.text('known: Z₂, Z₄, Z₆, …        unknown: the λᵢ themselves', W / 2, 420, {align: 'center', size: 19, alpha: seg(t, 11, 11.7)});
        K.ctx.restore();
      }
      // a 7 s pause at 31.3 explains the trial-state bound on Z24 while the ledger holds still
      // and a second one at 55 derives the purity bound before the ledger shows it
      const v = t < 31.3 ? t : t < 38.3 ? 31.3 : t < 55 ? t - 7 : t < 62 ? 48 : t - 14, tr = seg(t, 30.5, 31.2) * (1 - seg(t, 36.6, 37.2));
      if (b > 0) { K.ctx.save(); K.ctx.globalAlpha *= b; seedtoyBody(K, v - 14); K.ctx.restore(); }
      if (tr > 0) {
        const r = (s, y, at, o = {}) => K.text(s, 560, y, {size: 16, ...o, alpha: tr * seg(t, at, at + .6)});
        r('Z₂₄ = Σₘ e^(−βEₘ)  ≥  e^(−βE₀)', 66, 30.6);
        r('trial state ψ:   E₀ ≤ ⟨ψ|H|ψ⟩', 98, 31.8);
        r('(cheap for an MPS ψ, even on 24 sites)', 124, 32.4, {color: th.mute, size: 14});
        r('⟨ψ|H|ψ⟩ < 0  ⇒  e^(−βE₀) > 1  ⇒  Z₂₄ ≥ 1', 158, 33.5, {color: th.space});
        r('energies shifted by a constant per site; it cancels in the purities', 186, 34.5, {color: th.mute, size: 13});
      }
      const pd = seg(t, 55.3, 56) * (1 - seg(t, 68.4, 69));
      if (pd > 0) {
        K.tex('P_{\\mathrm s}(24)=\\frac{Z_{48}}{Z_{24}^{\\,2}}=\\frac{\\sum\\lambda^{48}}{\\big(\\sum\\lambda^{24}\\big)^2}', 740, 62, {align: 'center', size: 18, alpha: pd});
        K.tex('\\ge\\frac{\\lambda_0^{48}}{(\\lambda_0^{24}+0.0030)^2}', 740, 128, {align: 'center', size: 18, alpha: pd * seg(t, 57.5, 58.2)});
        K.text('top: keep λ₀ only;  bottom: rest ≤ 0.0030', 740, 172, {align: 'center', color: th.mute, size: 13, alpha: pd * seg(t, 58.5, 59.2)});
        K.tex('=\\big(1+0.0030/\\lambda_0^{24}\\big)^{-2}\\ \\ge\\ \\big(1+0.0030/0.97^{24}\\big)^{-2}', 740, 204, {align: 'center', size: 16, color: th.space, alpha: pd * seg(t, 62.5, 63.2)});
      }
    }
  });

  // The real seed, drawn as the closed tensor networks it evaluates: rings of X (plain and with a
  // symmetry twist g) for the moments, and the double-layer ring of a matrix-product trial state.
  const DBLV = [[36, .0479158, .1815202], [72, .0604996, .1733812], [144, .0686606, .1632381], [288, .0713616, .1447092],
    [576, .0632939, .1055161], [1152, .0375793, .0445941], [2304, .0084513, .0054986]];
  // a row of n boxes joined by legs, with the wrap-around leg drawn above (up) or below
  function boxRing(K, x0, dx, n, y, lab, o) {
    const {th} = K, x1 = x0 + (n - 1) * dx, h = o.up ? -34 : 34;
    K.path([[x1 + 30, y], [x1 + 45, y], [x1 + 45, y + h], [x0 - 45, y + h], [x0 - 45, y], [x0 - 30, y]], {color: th.ink, w: 2, alpha: o.alpha});
    K.line(x0 - 30, y, x1 + 30, y, {color: th.ink, w: 2, alpha: o.alpha});
    for (let i = 0; i < n; i++) K.box(x0 + i * dx, y, 36, 28, lab, {color: o.color, fill: o.fill, alpha: o.alpha});
  }
  chapter({
    id: 'seed', title: 'Where the starting numbers come from', dur: 54,
    beats: [[0, 'The real certificate (the inputs of OpenAI’s 2026 proof) runs the toy’s three steps on the spin-1 chain. Step 1 needs exact moments Zₙ = Tr Xⁿ: the closed ring of n boxes X, which is Tr e^(−βH) on a ring of n ≤ 12 sites, evaluated in exact integer arithmetic.'],
            [8, 'A sharper version inserts a symmetry operation g into the ring before closing it: Zₙ(g) = Tr(U_g Xⁿ). These twisted traces separate the eigenvalues by symmetry sector, so each sector gets its own moment bound.'],
            [16, 'Step 2 needs a trial state on 72 sites: a matrix-product state. One box A per site, with a physical leg (the spin) and two virtual legs of dimension 26 joined around the ring.'],
            [23, 'Its norm ⟨ψ|ψ⟩ glues ψ to its complex conjugate (the bra) along the physical legs: a double-layer ring, Tr T⁷² with T = Σₛ Aₛ ⊗ Āₛ, which is Σₛ Aₛ ⊗ Aₛ here because this MPS is real. Its energy ⟨ψ|H|ψ⟩ is 72 times the same ring with one h inserted. Both are rings of small matrices, contracted exactly.'],
            [31, 'Why that gives Z₇₂ > 1: energies carry a constant shift κ per site, which cancels in every purity, so Z₇₂ = Σₘ e^(−β(Eₘ + 72κ)) ≥ e^(−β(E₀ + 72κ)). The trial state has energy below −72κ, so E₀ < −72κ and that one factor already exceeds 1.'],
            [39, 'Then, as in the toy: Z₇₂ = Σ λ⁷² > 1 and the moment caps on every λ together force one eigenvalue λ₀ > 0.999. That certifies the seed at n₀ = 36: Pₛ at length 36 and Pτ at length 72.'],
            [45, 'The doubling step, run six times, carries them to (2304, 784), both defects below 0.01.']],
    draw(K, t) {
      const {th} = K;
      // step 1: rings of X, plain and twisted
      const P1 = seg(t, .3, 1) * (1 - seg(t, 15.4, 16));
      if (P1 > 0) {
        boxRing(K, 230, 90, 6, 170, 'X', {alpha: P1});
        K.box(455, 204, 28, 22, 'g', {color: th.time, fill: th.timeF, size: 15, alpha: P1 * seg(t, 8.5, 9.2)});
        K.tex('Z_n=\\operatorname{Tr}X^n=\\operatorname{Tr}e^{-\\beta H}\\ \\ (\\text{ring of }n\\le 12\\text{ sites})', W / 2, 280, {align: 'center', size: 20, alpha: P1 * seg(t, .8, 1.6)});
        K.text('evaluated exactly, in integers', W / 2, 322, {align: 'center', color: th.ok, size: 16, alpha: P1 * seg(t, 3, 3.7)});
        K.tex('Z_n(g)=\\operatorname{Tr}\\big(U_g\\,X^n\\big)', W / 2, 385, {align: 'center', size: 20, color: th.time, alpha: P1 * seg(t, 9.5, 10.2)});
        K.text('twists separate the eigenvalues by symmetry sector: one moment bound per sector', W / 2, 430, {align: 'center', color: th.mute, size: 16, alpha: P1 * seg(t, 11, 11.7)});
      }
      // step 2: the matrix-product trial state, its norm and its energy
      const P2 = seg(t, 16, 16.7) * (1 - seg(t, 44.4, 45)), x0 = 160, dx = 82, n = 8, yA = 130, yB = 250;
      if (P2 > 0) {
        const dbl = seg(t, 23.5, 24.3), en = seg(t, 27, 27.7);
        for (let i = 0; i < n; i++) {
          const x = x0 + i * dx, hcol = en > 0 && (i === 3 || i === 4);
          K.line(x, yA + 14, x, dbl > 0 ? yB - 14 : yA + 46, {color: th.space, w: 2, alpha: P2});
          if (dbl === 0) K.text('s', x + 8, yA + 50, {color: th.space, size: 13, italic: true, alpha: P2 * seg(t, 17, 17.7)});
          if (hcol && i === 3) K.box(x + dx / 2, (yA + yB) / 2, dx + 44, 30, 'h', {color: th.fail, fill: th.spaceF, alpha: P2 * en});
        }
        boxRing(K, x0, dx, n, yA, 'A', {up: true, alpha: P2});
        boxRing(K, x0, dx, n, yB, 'Ā', {alpha: P2 * dbl});
        K.text('virtual legs: dimension 26', x0 + n * dx + 20, yA - 6, {color: th.mute, size: 14, alpha: P2 * seg(t, 17.5, 18.2) * (1 - dbl)});
        K.text('physical leg: the spin', x0 + n * dx + 20, yA + 50, {color: th.space, size: 14, alpha: P2 * seg(t, 18.5, 19.2) * (1 - dbl)});
        K.text('ψ', x0 - 70, yA + 6, {size: 20, italic: true, alpha: P2});
        K.text('ψ̄', x0 - 70, yB + 6, {size: 20, italic: true, alpha: P2 * dbl});
        K.rect(x0 + dx - 26, yA - 22, 52, yB - yA + 44, {stroke: th.proved, w: 2, dash: [5, 4], alpha: P2 * seg(t, 25, 25.7) * (1 - en)});
        K.text('T', x0 + dx, yB + 42, {align: 'center', color: th.proved, size: 17, alpha: P2 * seg(t, 25, 25.7) * (1 - en)});
        K.tex('\\langle\\psi|\\psi\\rangle=\\operatorname{Tr}T^{72},\\qquad T=\\textstyle\\sum_s A_s\\otimes\\bar A_s', W / 2, 330, {align: 'center', size: 19, alpha: P2 * seg(t, 25.5, 26.2)});
        K.tex('\\langle\\psi|H|\\psi\\rangle=72\\times(\\text{ring with one }h)', W / 2, 372, {align: 'center', size: 19, color: th.fail, alpha: P2 * seg(t, 27.5, 28.2)});
        const w = seg(t, 31.5, 32.2);
        K.tex('Z_{72}=\\sum_m e^{-\\beta(E_m+72\\kappa)}\\ \\ge\\ e^{-\\beta(E_0+72\\kappa)}', W / 2, 420, {align: 'center', size: 19, alpha: P2 * w});
        K.tex('E_0\\le\\frac{\\langle\\psi|H|\\psi\\rangle}{\\langle\\psi|\\psi\\rangle}<-72\\kappa\\ \\Rightarrow\\ Z_{72}>1', W / 2, 474, {align: 'center', size: 19, color: th.space, alpha: P2 * seg(t, 34, 34.7)});
        K.text('+ moment bounds  ⇒  one λ₀ > 0.999  ⇒  seed n₀ = 36: Pₛ(36), Pτ(72)', W / 2, 520, {align: 'center', size: 17, color: th.ok, alpha: P2 * seg(t, 39.5, 40.2)});
      }
      // six doublings
      const B = seg(t, 45, 45.8);
      if (B > 0) {
        const gx = 120, gy = 70, gw = 700, gh = 340, lo = -2.4, hi = -.6, X = j => gx + j / 6 * gw, Y = v => gy + gh * (hi - Math.log10(v)) / (hi - lo);
        K.axes(gx, gy, gw, gh, {alpha: B});
        [[.12731, 'ε* ≈ 0.127', th.mute], [.01, 'target 0.01', th.proved]].forEach(([v, nm, col]) => { K.line(gx, Y(v), gx + gw, Y(v), {color: col, dash: [6, 5], alpha: B, w: 1.5}); K.text(nm, gx + gw + 8, Y(v) + 5, {color: col, size: 14, alpha: B}); });
        DBLV.forEach((r, j) => K.text(String(r[0]), X(j), gy + gh + 22, {align: 'center', color: th.mute, size: 13, alpha: B}));
        K.text('ring length (β doubles alongside, 12¼ → 784)', gx, gy + gh + 46, {color: th.mute, size: 14, alpha: B});
        const prog = ease(seg(t, 46, 52)) * 6;
        [[1, th.space, 'spatial defect'], [2, th.time, 'thermal defect']].forEach(([k, col, nm], q) => {
          const pts = DBLV.map((r, j) => [X(j), Y(r[k])]).filter((_, j) => j <= prog);
          K.path(pts, {color: col, w: 3, alpha: B}); pts.forEach(p => K.dot(p[0], p[1], 5, col, B));
          K.text(nm, gx + 12, gy + 22 + 22 * q, {color: col, size: 15, alpha: B});
        });
      }
    }
  });

  // Stretch a chapter's clock by S, so the densest figures play more slowly.
  function slow(list, id, S) {
    const c = list.find(x => x.id === id), d0 = c.draw;
    c.dur *= S; c.beats = c.beats.map(([b, txt]) => [b * S, txt]); c.draw = (K, t) => d0(K, t / S);
  }

  // Index of the caption beat showing at time t.
  function beatAt(ch, t) { let bi = 0; ch.beats.forEach((b, i) => { if (t >= b[0] - 1e-6) bi = i; }); return bi; }

  // Typeset TeX with MathJax into an image drawn at px font size; resolves to {img, w, h} in
  // logical pixels, rasterized at twice that resolution.
  function texImage(src, color, px) {
    return new Promise((resolve, reject) => {
      const svg = MathJax.tex2svg(src, {display: true}).querySelector('svg'), ex = px * .45;
      const w = parseFloat(svg.getAttribute('width')) * ex, h = parseFloat(svg.getAttribute('height')) * ex;
      svg.setAttribute('width', w * 2); svg.setAttribute('height', h * 2);
      const img = new Image(); img.onload = () => resolve({img, w, h}); img.onerror = reject;
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg.outerHTML.replace(/currentColor/g, color));
    });
  }

  root.STBV = {CH, kit, W, H, SERIF, clamp, seg, ease, lerp, slow, beatAt, texImage, minusH, realign, eigSym};
})(typeof window !== 'undefined' ? window : globalThis);
