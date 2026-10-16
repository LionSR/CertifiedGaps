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
  const SERIF = '"EB Garamond", Palatino, "DejaVu Serif", serif', SANS = SERIF, MONO = SERIF;

  // ---------- drawing kit ----------
  function kit(ctx, th, texImg) {
    const K = {ctx, th, W, H, seg, ease, lerp, clamp};
    K.text = (s, x, y, o = {}) => {
      ctx.save(); ctx.globalAlpha *= o.alpha ?? 1; ctx.fillStyle = o.color || th.ink;
      ctx.font = `${o.italic ? 'italic ' : ''}${o.weight || 400} ${Math.round((o.size || 16) * 1.12)}px ${o.font || SANS}`; ctx.textAlign = o.align || 'left';
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
      K.text('Δ', lx + 214, base - .45 * sc + 8, {size: 26, font: SERIF, color: th.proved, alpha: g});
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
      K.text('c', px + pw + 10, Y(.25) + 6, {size: 22, font: SERIF, color: th.proved, alpha: c});
      K.text('every spin-1 dot must stay above the dashed line', px, py + ph + 56, {color: th.proved, alpha: seg(t, 11.5, 12.5), size: 15});
    }
  });

  // 2 ------------------------------------------------------------------
  chapter({
    id: 'rate', title: 'A gap is a decay rate', dur: 14,
    beats: [[0, 'At inverse temperature β, level m has probability e^(−βEₘ)/Z.'],
            [4, 'Cooling down empties the excited levels. The ratio of the first two weights is exactly e^(−βΔ).'],
            [9.5, 'So it is enough to show: excited weight ≤ e^(−cβ) on every ring. Then Δ ≥ c.']],
    tex: ['\\frac{w_1}{w_0}=e^{-\\beta\\Delta}'],
    draw(K, t) {
      const {th} = K, E = [0, 1, 1, 1, 1.6, 2, 2.4], beta = lerp(.15, 4.5, ease(seg(t, 2, 9)));
      const w = E.map(e => Math.exp(-beta * e)), Z = w.reduce((x, y) => x + y);
      K.bars(70, 110, 330, 280, w.map(v => v / Z), th.time, {labels: true, alpha: seg(t, 0, .8)});
      K.text(`β = ${beta.toFixed(2)}`, 70, 440, {font: MONO, size: 18, alpha: seg(t, 0, .8)});
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
    tex: ['e^{-\\tau H}\\psi=\\sum_m c_m\\,e^{-\\tau E_m}|m\\rangle'],
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
        K.text(`E = ${e}`, x, y + R + 22, {align: 'center', color: th.mute, size: 13, font: MONO});
        const p = Math.pow(amp[m] * scale, 2), bh = 110;
        K.rect(x - 22, 400 - p * bh, 44, p * bh, {fill: cols[m], alpha: .85, r: 2});
        K.text(p.toFixed(2), x, 392 - p * bh, {align: 'center', color: th.mute, size: 12, font: MONO});
        if (t >= 18.5) K.text(`e^(−βE) = ${Math.exp(-3 * e).toFixed(3)}`, x, 432, {align: 'center', color: th.proved, size: 13, font: MONO, alpha: seg(t, 18.8, 19.5)});
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
  chapter({
    id: 'paths', title: 'A sum over histories', dur: 28,
    beats: [[0, 'Start with a fact about matrices. An entry of Aᴹ is a sum over paths through M steps: at each step pick an entry of A, multiply them all.'],
            [6, 'One spin that can flip, H = −σˣ. For a tiny step dτ, A = e^(−dτ H) ≈ [[1, dτ], [dτ, 1]]: staying costs 1, flipping costs dτ. And e^(−βH) = Aᴹ with M = β/dτ.'],
            [12.5, 'Make the steps tiny. A path becomes a world line that mostly stays and flips a few times. With k flips its weight is dτᵏ, and there are about (β/dτ)ᵏ/k! places for them: βᵏ/k! in total.'],
            [19, 'The trace keeps only paths that end where they started. So k is even, and the time axis closes into a circle.'],
            [23, 'Add all paths: Σ 2βᵏ/k! over even k = 2 cosh β, the exact answer. This is the Euclidean path integral: quantum mechanics as a sum over histories.']],
    tex: ['A=\\begin{pmatrix}1&d\\tau\\\\ d\\tau&1\\end{pmatrix}', 'Z=\\sum_{k\\ \\mathrm{even}}\\frac{2\\beta^k}{k!}=2\\cosh\\beta'],
    draw(K, t) {
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
        K.text('e^(−βH) = Aᴹ,  M = β / dτ', 700, 350, {color: th.ink, size: 15, font: MONO, alpha: A1 * seg(t, 9, 10)});
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
          K.text(`k = ${fl.length}:  weight βᵏ/k! = ${(Math.pow(2, fl.length) / [1, 1, 2, 6, 24][fl.length]).toFixed(2)} (β = 2)`, x1, yD + 26, {align: 'right', color: th.mute, size: 13, font: MONO, alpha: al});
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
        terms.forEach((v, i) => { if (i >= nshow) return; K.rect(X(i), Y(run + v), 50, Y(run) - Y(run + v), {fill: th.space, alpha: C * .85}); K.text(`k=${2 * i}`, X(i) + 25, gy + gh + 20, {align: 'center', color: th.mute, size: 12, font: MONO, alpha: C}); run += v; });
        K.line(gx, Y(2 * Math.cosh(beta)), gx + gw, Y(2 * Math.cosh(beta)), {color: th.proved, dash: [6, 4], w: 2, alpha: C});
        K.text(`2 cosh 2 = ${(2 * Math.cosh(2)).toFixed(3)}`, gx + gw, Y(2 * Math.cosh(beta)) - 8, {align: 'right', color: th.proved, size: 14, alpha: C});
        K.text(`running sum = ${run.toFixed(3)}`, gx, gy - 14, {font: MONO, size: 15, alpha: C});
        K.tex('Z=\\sum_{k\\ \\mathrm{even}}\\frac{2\\beta^k}{k!}=2\\cosh\\beta', 780, 260, {align: 'center', size: 24, alpha: C * seg(t, 25, 26)});
      }
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
  const CEV = [[0, [.15, .62]], [1, [.3, .8]], [2, [.12, .47, .7]], [3, [.25, .55]], [4, [.4, .88]], [5, [.2, .66]]];
  chapter({
    id: 'chain', title: 'Histories of a chain', dur: 26,
    beats: [[0, 'Warm-up: a classical Ising ring. Each site carries a spin ±1, and a 2×2 matrix T passes it to the next site. Then Z = Tr Tᴸ.'],
            [5, 'Quantum chain: every site gets a world line in imaginary time. Each bond term S·S fires events at random times; an event acts on both sites of its bond.'],
            [10.5, 'Now what is passed along the ring is the list of events on a bond: a on the left of a site, b on its right.'],
            [15, 'X(a, b): multiply that site’s matrices in time order and take the trace. With G = iS (real 3×3 matrices): no events gives 3; one z-event on each side gives Tr(G_z G_z) = −2.'],
            [21, 'So Z = Tr Xᴸ, exactly like Ising, only the index is a whole history instead of ±1. This X is the transfer matrix of the proof.']],
    draw(K, t) {
      const {th} = K, A = 1 - seg(t, 4.6, 5.2);
      if (A > 0) {
        const n = 6, x0 = 150, dx = 120, y = 220;
        for (let i = 0; i < n; i++) {
          const x = x0 + i * dx, up = [1, 1, -1, 1, -1, 1][i] > 0;
          K.circle(x, y, 26, {color: th.ink, w: 2, alpha: A}); K.arrow(x, y + (up ? 14 : -14), x, y + (up ? -14 : 14), {color: up ? th.space : th.time, w: 3, alpha: A, head: 8});
          K.rect(x + 36, y - 16, dx - 72, 32, {fill: th.provedF, stroke: th.proved, w: 1.5, r: 6, alpha: A * seg(t, 1.5, 2.2)});
          K.text('T', x + dx / 2, y + 6, {align: 'center', color: th.proved, size: 16, font: SERIF, alpha: A * seg(t, 1.5, 2.2)});
        }
        K.text('spins ±1 on sites; T(s, s′) = e^(K s s′) on bonds; the last bond wraps around', W / 2, 310, {align: 'center', color: th.mute, size: 15, alpha: A});
        K.text('Z = Tr Tᴸ', W / 2, 370, {align: 'center', size: 26, font: SERIF, alpha: A * seg(t, 2.5, 3.3)});
      }
      const B = seg(t, 5.2, 5.8);
      if (B > 0) {
        const N = 6, x0 = 110, dx = 88, y0 = 80, y1 = 400, Hh = y1 - y0, X = i => x0 + i * dx, hi = seg(t, 10.5, 11.3);
        K.rect(X(2) - 44, y0 - 6, 88, Hh + 12, {fill: th.spaceF, alpha: hi * B});
        for (let i = 0; i < N; i++) { K.line(X(i), y0, X(i), y1, {color: th.mute, w: 2, alpha: B}); K.text(`${i + 1}`, X(i), y1 + 22, {align: 'center', color: th.mute, size: 13, alpha: B}); }
        CEV.forEach(([bd, ts], k) => ts.forEach((f, q) => {
          const al = B * seg(t, 6 + (k * 2 + q) * .2, 6.5 + (k * 2 + q) * .2), yy = y1 - f * Hh;
          const col = hi > 0 && bd === 1 ? th.proved : hi > 0 && bd === 2 ? th.time : th.space;
          if (bd < N - 1) K.line(X(bd), yy, X(bd + 1), yy, {color: col, w: 3, alpha: al});
          else { K.line(X(bd), yy, X(bd) + 40, yy, {color: col, w: 3, alpha: al, dash: [5, 3]}); K.line(X(0) - 40, yy, X(0), yy, {color: col, w: 3, alpha: al, dash: [5, 3]}); }
          K.dot(X(bd), yy, 4, col, al); K.dot(X((bd + 1) % N), yy, 4, col, al);
        }));
        K.text('site', X(0) - 30, y1 + 22, {align: 'right', color: th.mute, size: 13, alpha: B});
        K.text('τ ↑', X(0) - 30, y0 + 10, {align: 'right', color: th.time, size: 14, alpha: B});
        K.text('a', X(1) + dx / 2, y0 - 12, {align: 'center', color: th.proved, size: 22, font: SERIF, alpha: hi * B});
        K.text('b', X(2) + dx / 2, y0 - 12, {align: 'center', color: th.time, size: 22, font: SERIF, alpha: hi * B});
        const c = seg(t, 15, 15.8);
        if (c > 0) {
          const px = 640, rows = [['a', 'b', 'X(a, b)'], ['none', 'none', 'Tr 1 = 3'], ['z', 'none', 'Tr G_z = 0'], ['z', 'z', 'Tr G_zG_z = −2'], ['z', 'x', 'Tr G_zG_x = 0']];
          rows.forEach((r, i) => { const yy = 120 + i * 40, al = c * seg(t, 15 + i * .5, 15.6 + i * .5);
            K.text(r[0], px, yy, {color: i ? th.proved : th.mute, size: 15, font: MONO, alpha: al}); K.text(r[1], px + 70, yy, {color: i ? th.time : th.mute, size: 15, font: MONO, alpha: al});
            K.text(r[2], px + 140, yy, {color: i ? th.ink : th.mute, size: 15, font: MONO, alpha: al});
            if (i === 0) K.line(px, yy + 10, px + 300, yy + 10, {color: th.line, w: 1, alpha: al}); });
          K.text('entries can be negative: that is fine', px, 340, {color: th.mute, size: 14, alpha: c * seg(t, 18, 19)});
        }
        K.text('Z = Tr Xᴸ', 760, 420, {align: 'center', size: 28, font: SERIF, color: th.space, alpha: seg(t, 21, 22)});
      }
    }
  });

  // 4 ------------------------------------------------------------------
  chapter({
    id: 'slice', title: 'Two ways to slice the torus', dur: 16,
    beats: [[0, 'Slice across: each row is a step of time evolution. Summing over energies gives Z = Σ e^(−βE).'],
            [5, 'Slice up: each column is one site carried along the ring by a transfer matrix X. Then Z = Tr Xᴸ = Σ λᴸ.'],
            [10.5, 'Same number, two sums. On a long ring the largest λ wins, just as at low temperature the lowest E wins.']],
    tex: ['Z=\\sum_m e^{-\\beta E_m}', 'Z=\\operatorname{Tr}X^L=\\sum_i\\lambda_i^{\\,L}'],
    draw(K, t) {
      const {th} = K, x = 70, y = 120, c = 40, cols = 8, rows = 6;
      K.torus(x, y, cols, rows, c);
      K.text('L sites', x + cols * c / 2, y + rows * c + 28, {align: 'center', color: th.space, size: 15});
      K.text('β', x - 18, y + rows * c / 2 + 6, {align: 'center', color: th.time, size: 18, font: SERIF});
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
    beats: [[0, 'A state with probabilities p₁, p₂, … has purity P = Σ p². Draw each p as a square: P is their total area.'],
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
      K.text(`P = Σp² = ${P.toFixed(3)}`, 600, 200, {font: MONO, size: 22, alpha: seg(t, 2, 3)});
      K.text(`p₀ = ${p[0].toFixed(3)}`, 600, 240, {font: MONO, size: 22, alpha: r});
      K.text('squares ⊂ dashed box  ⇒  P ≤ p₀', 600, 300, {size: 18, color: th.proved, alpha: seg(t, 6.5, 7.5)});
    }
  });

  // 6 ------------------------------------------------------------------
  chapter({
    id: 'two', title: 'Two purities from one torus', dur: 13,
    beats: [[0, 'Normalize the torus in space or in time and you get two density matrices, σ and ρ.'],
            [3, 'The purity of σ glues two copies side by side: the spatial purity Pₛ.'],
            [7.5, 'The purity of ρ stacks two copies: the thermal purity Pτ. Both are ratios of partition functions a computer can evaluate.']],
    tex: ['P_{\\mathrm s}(n,\\beta)=\\frac{Z_{2n}(\\beta)}{Z_n(\\beta)^2}', 'P_\\tau(n,\\beta)=\\frac{Z_n(2\\beta)}{Z_n(\\beta)^2}'],
    draw(K, t) {
      const {th} = K, c = 26;
      const m1 = ease(seg(t, 3.5, 5.5)), x1 = 90, y1 = 200;
      K.torus(x1, y1, 5, 3, c, {fill: th.spaceF}); K.torus(x1 + 5 * c + lerp(60, 0, m1), y1, 5, 3, c, {fill: th.spaceF});
      K.text('side by side', x1, y1 - 20, {color: th.space, size: 16});
      K.tex('P_{\\mathrm s}(n,\\beta)=\\frac{Z_{2n}(\\beta)}{Z_n(\\beta)^2}', x1 + 140, 390, {align: 'center', size: 24, color: th.space, alpha: seg(t, 5, 6)});
      const m2 = ease(seg(t, 8, 10)), x2 = 640, y2 = 110;
      K.torus(x2, y2, 5, 3, c, {fill: th.timeF}); K.torus(x2, y2 + 3 * c + lerp(50, 0, m2), 5, 3, c, {fill: th.timeF});
      K.text('stacked', x2, y2 - 20, {color: th.time, size: 16});
      K.tex('P_\\tau(n,\\beta)=\\frac{Z_n(2\\beta)}{Z_n(\\beta)^2}', x2 + 65, 390, {align: 'center', size: 24, color: th.time, alpha: seg(t, 9.5, 10.5)});
    }
  });

  // 7 ------------------------------------------------------------------
  const fsq = z => z * z / (2 * (1 - z) * (1 - z));
  chapter({
    id: 'square', title: 'Squaring sharpens', dur: 15,
    beats: [[0, 'Square every probability and renormalize. The leader pulls ahead.'],
            [5, 'Gluing two copies of a torus is exactly this squaring. A defect ε = 1 − P becomes about ε²/2.'],
            [10, 'Start below ½ and repeated squaring drives the defect to zero: the cobweb drops off the diagonal.']],
    tex: ['f(z)=\\frac{z^2}{2(1-z)^2}'],
    draw(K, t) {
      const {th} = K;
      let p = [.5, .3, .2]; const k = seg(t, .8, 8) * 3, n = Math.floor(k), fr = ease(k - n);
      const sq = q => { const s = q.reduce((u, v) => u + v * v, 0); return q.map(v => v * v / s); };
      for (let i = 0; i < Math.min(n, 3); i++) p = sq(p);
      const nx = n < 3 ? sq(p) : p, cur = p.map((v, i) => lerp(v, nx[i], n < 3 ? fr : 0));
      K.bars(70, 110, 300, 300, cur, th.ink, {labels: true});
      const P = cur.reduce((u, v) => u + v * v, 0);
      K.text(`squarings: ${Math.min(n, 3)}`, 70, 450, {font: MONO, size: 16});
      K.text(`defect 1 − P = ${(1 - P).toFixed(4)}`, 70, 478, {font: MONO, size: 16, color: th.proved});
      const gx = 480, gy = 100, gw = 410, gh = 330, zm = .5, X = z => gx + z / zm * gw, Y = z => gy + gh - z / zm * gh, b = seg(t, 5, 5.8);
      K.axes(gx, gy, gw, gh, {alpha: b});
      K.line(X(0), Y(0), X(zm), Y(zm), {color: th.mute, dash: [4, 4], alpha: b, w: 1.5});
      const fp = []; for (let z = 0; z <= zm + 1e-9; z += .005) fp.push([X(z), Y(fsq(z))]); K.path(fp, {color: th.proved, w: 3, alpha: b});
      K.tex('f(z)=\\frac{z^2}{2(1-z)^2}', gx + 30, gy + 30, {size: 20, color: th.proved, alpha: b});
      K.text('old defect →', gx + gw, gy + gh + 24, {align: 'right', color: th.mute, size: 13, alpha: b});
      K.text('new defect', gx - 8, gy - 8, {color: th.mute, size: 13, alpha: b});
      let z = .45; const pts = [[X(z), Y(0)]]; for (let i = 0; i < 5; i++) { const fz = fsq(z); pts.push([X(z), Y(fz)], [X(fz), Y(fz)]); z = fz; }
      const cw = seg(t, 10, 14.5) * (pts.length - 1), ni = Math.floor(cw), sh = pts.slice(0, ni + 1);
      if (ni < pts.length - 1) sh.push([lerp(pts[ni][0], pts[ni + 1][0], cw - ni), lerp(pts[ni][1], pts[ni + 1][1], cw - ni)]);
      K.path(sh, {color: th.ink, w: 2});
    }
  });

  // 8 ------------------------------------------------------------------
  chapter({
    id: 'couple', title: 'Space and time are tied', dur: 15,
    beats: [[0, 'Squaring grows the torus in one direction at a time. We need both.'],
            [2.5, 'Glue four copies into one twice as long and twice as tall: rows first …'],
            [6.5, '… or columns first. Both give the same torus.'],
            [10.5, 'Writing each step as a purity gives one identity linking the spatial and thermal purities.']],
    tex: ['P_{\\mathrm s}(n,2\\beta)\\,P_\\tau(n,\\beta)^2=P_{\\mathrm s}(n,\\beta)^2\\,P_\\tau(2n,\\beta)'],
    draw(K, t) {
      const {th} = K, c = 14, L = 150, R = 600, T = 70, B = 290;
      const g1 = ease(seg(t, 3, 4.5)), g2 = ease(seg(t, 4.5, 6)), h1 = ease(seg(t, 7, 8.5)), h2 = ease(seg(t, 8.5, 10));
      const tiles = (x, y, cols, rows, nx, ny, al, fill) => { for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) K.torus(x + i * (cols * c + 14), y + j * (rows * c + 14), cols, rows, c, {alpha: al, fill}); };
      tiles(L, T, 4, 3, 2, 2, 1);
      tiles(R, T, 8, 3, 1, 2, Math.max(.2, g1));
      tiles(L, B, 4, 6, 2, 1, Math.max(.2, h1));
      K.torus(R, B, 8, 6, c, {alpha: Math.max(.2, g2, h2), fill: th.provedF});
      K.arrow(L + 140, T + 30, R - 20, T + 30, {color: th.space, w: 3, alpha: lerp(.25, 1, g1)});
      K.arrow(R + 56, T + 104, R + 56, B - 12, {color: th.time, w: 3, alpha: lerp(.25, 1, g2)});
      K.arrow(L + 56, T + 104, L + 56, B - 12, {color: th.time, w: 3, alpha: lerp(.25, 1, h1)});
      K.arrow(L + 140, B + 42, R - 20, B + 42, {color: th.space, w: 3, alpha: lerp(.25, 1, h2)});
      K.text('× Pₛ(n,β)²', (L + R) / 2 + 70, T + 20, {align: 'center', color: th.space, size: 15, alpha: lerp(.25, 1, g1)});
      K.text('× Pτ(2n,β)', R + 70, (T + B) / 2 + 50, {color: th.time, size: 15, alpha: lerp(.25, 1, g2)});
      K.text('× Pτ(n,β)²', L + 42, (T + B) / 2 + 50, {align: 'right', color: th.time, size: 15, alpha: lerp(.25, 1, h1)});
      K.text('× Pₛ(n,2β)', (L + R) / 2 + 70, B + 32, {align: 'center', color: th.space, size: 15, alpha: lerp(.25, 1, h2)});
      K.tex('P_{\\mathrm s}(n,2\\beta)\\,P_\\tau(n,\\beta)^2=P_{\\mathrm s}(n,\\beta)^2\\,P_\\tau(2n,\\beta)', W / 2, 490, {align: 'center', size: 22, alpha: seg(t, 11, 12)});
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
            [14, 'Every purity lies between 0 and 1 and climbs to 1. This is the behaviour the bootstrap relies on.'],
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
        K.text('pᵢ = λᵢⁿ / Zₙ', W / 2, 230, {align: 'center', size: 34, font: SERIF, alpha: seg(t, .5, 1.5)});
        K.text('must be probabilities: real, ≥ 0, adding to 1', W / 2, 280, {align: 'center', size: 18, italic: true, alpha: seg(t, 2, 3)});
      }
    }
  });

  chapter({
    id: 'rpwhere', title: 'Where the symmetry comes from', dur: 38,
    beats: [[0, 'Recall X(a, b): one site’s matrices multiplied in time order, then traced. a lists the events that arrive from the left bond, b those from the right bond.'],
            [7, 'Heisenberg: an event puts the same real matrix G on the site whichever side it comes from. Example: a = {G_z at 0.2, G_x at 0.6}, b = {G_y at 0.4}.'],
            [14, 'In time order the site sees G_z, G_y, G_x, so X(a, b) = Tr(G_x G_y G_z) = 1. Swap the two sides: the site sees the same three matrices at the same times, so X(b, a) = 1 too.'],
            [22, 'Counter-example: a bond that puts G_z on its left site and the projector P = diag(1, 0, 0) on its right site. Now the matrix depends on the side.'],
            [28, 'Two events from the left give Tr(P P) = 1; the same two events from the right give Tr(G_z G_z) = −2. Not symmetric. Reflection positivity is the condition that rules this out: the bond must look the same from both sides, with positive weights.']],
    draw(K, t) {
      const {th} = K, sx = 300, y0 = 110, y1 = 400, Yt = f => y1 - f * (y1 - y0);
      const A = 1 - seg(t, 21.6, 22.2);
      if (A > 0) {
        const sw = ease(seg(t, 17.5, 19.5));
        K.line(sx, y0, sx, y1, {color: th.ink, w: 3, alpha: A});
        K.text('one site', sx, y1 + 26, {align: 'center', color: th.mute, size: 14, alpha: A});
        K.text('from the left bond', sx - 160, y0 - 16, {align: 'center', color: th.proved, size: 14, alpha: A});
        K.text('from the right bond', sx + 160, y0 - 16, {align: 'center', color: th.time, size: 14, alpha: A});
        K.text('τ ↑', sx - 260, y0 + 6, {color: th.time, size: 14, alpha: A});
        const ev = [[.2, -1, 'z'], [.6, -1, 'x'], [.4, 1, 'y']];
        ev.forEach(([f, s0, g], i) => {
          const s = lerp(s0, -s0, sw), y = Yt(f), al = A * seg(t, 7.5 + i * .8, 8.2 + i * .8), col = s0 < 0 ? th.proved : th.time;
          K.line(sx + s * 140, y, sx, y, {color: col, w: 3, alpha: al}); K.dot(sx + s * 140, y, 5, col, al); K.dot(sx, y, 5, th.ink, al);
          K.text(`G_${g}`, sx + (s > 0 ? 14 : -14), y - 8, {align: s > 0 ? 'left' : 'right', size: 16, alpha: al});
          K.text(`τ = ${f}`, sx + s * 150 + (s > 0 ? 6 : -6), y + 5, {align: s > 0 ? 'left' : 'right', color: th.mute, size: 13, alpha: al});
        });
        K.text('time-ordered product (later on the left):', 560, 170, {size: 16, alpha: A * seg(t, 14, 14.8)});
        K.text('X = Tr(G_x G_y G_z) = 1', 560, 205, {size: 22, font: SERIF, alpha: A * seg(t, 15, 15.8)});
        K.text(sw > .5 ? 'sides swapped: same matrices, same times' : '', 560, 260, {size: 16, color: th.ok, alpha: A * seg(t, 19.5, 20.2)});
        K.text('X(a, b) = X(b, a) = 1', 560, 295, {size: 22, font: SERIF, color: th.ok, alpha: A * seg(t, 20, 20.8)});
      }
      const B = seg(t, 22.2, 23);
      if (B > 0) {
        K.text('a bond that is not the same on both sides:', 80, 90, {size: 17, alpha: B});
        K.text('G_z on the left site,  P = diag(1, 0, 0) on the right site', 80, 120, {size: 17, color: th.fail, alpha: B});
        const panel = (x0, side, label, val, at) => {
          const al = B * seg(t, at, at + .8), x = x0 + 120;
          K.line(x, 170, x, 400, {color: th.ink, w: 3, alpha: al});
          [.3, .7].forEach(f => { const y = 400 - f * 230; K.line(x + side * 120, y, x, y, {color: side < 0 ? th.proved : th.time, w: 3, alpha: al}); K.dot(x, y, 5, th.ink, al);
            K.text(label, x + (side > 0 ? 14 : -14), y - 8, {align: side > 0 ? 'left' : 'right', size: 16, alpha: al}); });
          K.text(side < 0 ? 'two events from the left' : 'the same events from the right', x, 430, {align: 'center', color: th.mute, size: 14, alpha: al});
          K.text(val, x, 470, {align: 'center', size: 22, font: SERIF, alpha: B * seg(t, at + 3, at + 3.8)});
        };
        panel(140, -1, 'P', 'Tr(P P) = 1', 23);
        panel(520, 1, 'G_z', 'Tr(G_z G_z) = −2', 25);
        K.text('X(a, ∅) = 1  ≠  X(∅, a) = −2', W / 2, 515, {align: 'center', size: 18, color: th.fail, alpha: seg(t, 30, 31)});
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
  chapter({
    id: 'rptest', title: 'Testing a bond for reflection positivity', dur: 28,
    beats: [[0, 'A general test. Write −h for two spin-1 sites as a 9 × 9 matrix with entries ⟨i j| −h |k l⟩. Here i, k belong to the left site (gold) and j, l to the right site (blue).'],
            [5.5, 'Follow one entry: ⟨x x| −h |y y⟩ = 1. Move it to row (i, k) = (x, y), column (j, l) = (x, y). Do this for every entry: the matrix is “realigned”.'],
            [12, 'Why: realignment turns a product A ⊗ B into an outer product vec(A) vec(B)ᵀ. So −h = Σ G_α ⊗ G_α becomes a sum of squares Σ v vᵀ, with v = vec(G_α).'],
            [18, 'A sum of squares has eigenvalues ≥ 0. Here they are 2, 2, 2, 0, …: the antiferromagnet passes.'],
            [22.5, 'Ferromagnet: −h = −Σ G_α ⊗ G_α, eigenvalues −2. Not a sum of squares, not reflection positive. In general: realigned −h ⪰ 0, up to an on-site shift, means reflection positive.']],
    draw(K, t) {
      const {th} = K, fm = t >= 22.5, M = minusH(fm ? -1 : 1), cs = 30, ox = 140, oy = 80;
      const mv = ease(seg(t, 8.5, 11.5)), one = ease(seg(t, 6.5, 8.5));
      const labelsH = r => [Math.floor(r / 3), r % 3];
      // row/column labels
      for (let r = 0; r < 9; r++) {
        const [p, q] = labelsH(r), al = 1;
        const rowL = mv < .5 ? [IX[p], th.proved, IX[q], th.time] : [IX[p], th.proved, IX[q], th.proved];
        const colL = mv < .5 ? [IX[p], th.proved, IX[q], th.time] : [IX[p], th.time, IX[q], th.time];
        K.text(rowL[0], ox - 26, oy + r * cs + cs / 2 + 5, {align: 'center', color: rowL[1], size: 12, font: MONO, alpha: al});
        K.text(rowL[2], ox - 14, oy + r * cs + cs / 2 + 5, {align: 'center', color: rowL[3], size: 12, font: MONO, alpha: al});
        K.text(colL[0], ox + r * cs + cs / 2 - 5, oy - 8, {align: 'center', color: colL[1], size: 12, font: MONO, alpha: al});
        K.text(colL[2], ox + r * cs + cs / 2 + 5, oy - 8, {align: 'center', color: colL[3], size: 12, font: MONO, alpha: al});
      }
      K.text(mv < .5 ? 'rows (i, j)   columns (k, l)' : 'rows (i, k)   columns (j, l)', ox, oy + 9 * cs + 28, {color: th.mute, size: 14, font: MONO});
      K.rect(ox - 2, oy - 2, 9 * cs + 4, 9 * cs + 4, {stroke: th.line, w: 1.5});
      for (let a = 1; a < 3; a++) { K.line(ox + a * 3 * cs, oy, ox + a * 3 * cs, oy + 9 * cs, {color: th.line, w: 1}); K.line(ox, oy + a * 3 * cs, ox + 9 * cs, oy + a * 3 * cs, {color: th.line, w: 1}); }
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) for (let l = 0; l < 3; l++) {
        const v = M[i * 3 + j][k * 3 + l]; if (Math.abs(v) < 1e-9) continue;
        const special = i === 0 && j === 0 && k === 1 && l === 1;
        const f = special ? Math.max(one, mv) : mv;
        const r = lerp(i * 3 + j, i * 3 + k, f), c = lerp(k * 3 + l, j * 3 + l, f);
        K.rect(ox + c * cs + 2, oy + r * cs + 2, cs - 4, cs - 4, {fill: v > 0 ? th.proved : th.time, r: 4, alpha: special && t >= 5.5 && t < 12 ? 1 : .8});
        if (special && t >= 5.5 && t < 12) K.rect(ox + c * cs, oy + r * cs, cs, cs, {stroke: th.ink, w: 2.5, r: 5});
        K.text(v > 0 ? '+' : '−', ox + c * cs + cs / 2, oy + r * cs + cs / 2 + 5, {align: 'center', size: 13, font: MONO, color: th.bg});
      }
      if (t >= 5.5 && t < 12) K.text('⟨x x| −h |y y⟩ = +1  →  row (x, y), column (x, y)', 470, 120, {size: 15, font: MONO, alpha: seg(t, 5.8, 6.5)});
      const D = seg(t, 12, 12.8) * (1 - seg(t, 17.6, 18.2));
      if (D > 0) {
        K.text('realigned −h  =', 470, 150, {size: 17, alpha: D});
        for (let a = 0; a < 3; a++) {
          const v = []; for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) v.push(G3[a][i][k]);
          const bx = 470 + a * 150, by = 180, c2 = 13;
          for (let r = 0; r < 9; r++) for (let cc = 0; cc < 9; cc++) { const x = v[r] * v[cc] * (fm ? -1 : 1); if (!x) continue; K.rect(bx + cc * c2, by + r * c2, c2 - 1, c2 - 1, {fill: x > 0 ? th.proved : th.time, alpha: D * seg(t, 12.5 + a * .6, 13.1 + a * .6)}); }
          K.rect(bx - 1, by - 1, 9 * c2 + 1, 9 * c2 + 1, {stroke: th.line, w: 1, alpha: D});
          K.text(`vec(G_${IX[a]}) vec(G_${IX[a]})ᵀ`, bx + 58, by + 9 * c2 + 22, {align: 'center', size: 12, font: MONO, color: th.mute, alpha: D});
          if (a < 2) K.text('+', bx + 9 * c2 + 15, by + 64, {align: 'center', size: 22, alpha: D});
        }
        K.text('each term is a square (v vᵀ): never negative', 470, 380, {color: th.ok, size: 15, alpha: D * seg(t, 14.5, 15.3)});
      }
      const E2 = seg(t, 18, 18.8);
      if (E2 > 0) {
        const ev = eigSym(realign(M)), gx = 500, gy = 140, gw = 380, gh = 230, Y = v => gy + gh / 2 - v / 2.4 * gh / 2, bw = gw / 9 - 6;
        K.line(gx, Y(0), gx + gw, Y(0), {color: th.mute, w: 1.5, alpha: E2});
        ev.forEach((v, i) => K.rect(gx + i * (bw + 6), Math.min(Y(0), Y(v)), bw, Math.max(2, Math.abs(Y(v) - Y(0))), {fill: v < -1e-9 ? th.fail : th.ok, alpha: E2, r: 2}));
        K.text('eigenvalues of the realigned matrix', gx, gy - 14, {color: th.mute, size: 14, alpha: E2});
        K.text(fm ? 'ferromagnet: −2, −2, −2  ✗ not reflection positive' : 'antiferromagnet: 2, 2, 2, 0, …  ✓ reflection positive', gx, gy + gh + 34, {color: fm ? th.fail : th.ok, size: 16, alpha: E2});
      }
    }
  });

  // 10 -----------------------------------------------------------------
  chapter({
    id: 'double', title: 'One doubling step', dur: 17,
    beats: [[0, 'Start with both purities at least 1 − ε at size (n, β).'],
            [2.5, 'Move 1: the coupling identity turns them into a spatial bound at 2β, at the cost of tripling ε.'],
            [6, 'Move 2: square in space. The defect drops to about 5ε².'],
            [9.5, 'Moves 3 and 4: the same in the thermal lane, using move 2’s result.'],
            [13.5, 'Out comes the same pair of bounds at (2n, 2β) with ε replaced by 5ε². The step repeats.']],
    draw(K, t) {
      const {th} = K, e = .01, d = .5 * Math.pow(Math.pow(1 - e, -3) - 1, 2);
      const xs = 90, xt = 520, bw = 350, bh = 50, rows = [70, 160, 250, 340, 430];
      const box = (x, r, col, s, at, dash) => { const al = seg(t, at, at + .7); K.rect(x, rows[r], bw, bh, {fill: th.surface, stroke: col, w: 2, r: 8, dash, alpha: al}); K.text(s, x + 16, rows[r] + 32, {font: MONO, size: 15, alpha: al}); };
      box(xs, 0, th.space, 'Pₛ(n,β) ≥ 1 − ε = 0.99', .3); box(xt, 0, th.time, 'Pτ(2n,β) ≥ 1 − ε = 0.99', .3);
      box(xs, 1, th.space, `Pₛ(n,2β) ≥ (1−ε)³ = ${Math.pow(1 - e, 3).toFixed(4)}`, 3, [6, 5]);
      box(xs, 2, th.space, `Pₛ(2n,2β) ≥ 1 − ${d.toFixed(5)}`, 6.5);
      box(xt, 3, th.time, `Pτ(4n,β) ≥ ${(Math.pow(1 - e, 2) * (1 - d)).toFixed(4)}`, 10, [6, 5]);
      box(xt, 4, th.time, 'Pτ(4n,2β) ≥ 1 − 5ε² = 0.9995', 12);
      const mv = (x, r1, r2, k, lab, at) => { const al = seg(t, at - .8, at); K.arrow(x, rows[r1] + bh, x, rows[r2] - 2, {color: th.proved, w: 2.5, alpha: al}); K.text(`${k}  ${lab}`, x + 12, (rows[r1] + bh + rows[r2]) / 2 + 5, {color: th.proved, size: 14, alpha: al}); };
      mv(xs + 40, 0, 1, '①', 'couple', 3); mv(xs + 40, 1, 2, '②', 'square', 6.5); mv(xt + 40, 0, 3, '③', 'couple', 10); mv(xt + 40, 3, 4, '④', 'square', 12);
      K.arrow(xs + bw, rows[2] + bh / 2, xt - 6, rows[3] + 10, {color: th.mute, w: 1.5, alpha: seg(t, 9.5, 10)});
      K.arrow(xt, rows[0] + bh / 2, xs + bw - 30, rows[1] + 4, {color: th.mute, w: 1.5, alpha: seg(t, 2.5, 3)});
      K.text('(n, β, ε)  →  (2n, 2β, 5ε²)', xs, 515, {font: SERIF, size: 24, color: th.proved, alpha: seg(t, 13.5, 14.5)});
    }
  });

  // 11 -----------------------------------------------------------------
  chapter({
    id: 'iterate', title: 'Repeat, and read off the gap', dur: 16,
    beats: [[0, 'Each step doubles β and roughly squares ε. Plot log ε against β.'],
            [6, 'The points fall on a straight line: ε decays exponentially in β.'],
            [10.5, 'An exponential decay rate in β is exactly a gap. For spin 1: Δ(L) > log 20 / 784 ≈ 0.0038 for every even L ≥ 2304.']],
    draw(K, t) {
      const {th} = K, gx = 110, gy = 60, gw = 760, gh = 360, bmax = 13000, lmin = -24;
      const X = b => gx + b / bmax * gw, Y = l => gy + l / lmin * gh;
      K.axes(gx, gy, gw, gh);
      [0, -5, -10, -15, -20].forEach(l => { K.line(gx, Y(l), gx + gw, Y(l), {color: th.line, w: 1}); K.text(l ? `10${String(l).split('').map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(c)]).join('')}` : '1', gx - 8, Y(l) + 5, {align: 'right', color: th.mute, size: 13, font: MONO}); });
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
      K.text('Δ(L) > log 20 / 784', X(8200), Y(-5), {color: th.proved, size: 24, font: SERIF, alpha: seg(t, 11, 12)});
    }
  });

  // the seed certificate, run on a toy spectrum with every number shown
  chapter({
    id: 'seedtoy', title: 'The starting numbers, worked on a toy', dur: 46,
    beats: [[0, 'The certificate never sees the eigenvalues of X. To see how it still bounds them, hide four: 1, 0.6, −0.5, 0.4. A computer knows only exact traces Zₙ of small rings.'],
            [8, 'Step 1, moments. Z₆ = Σλ⁶ = 1.0664, computed exactly. Every term λ⁶ is ≥ 0, so no single |λ| can exceed 1.0664^(1/6) = 1.011.'],
            [16, 'Step 2, a trial state. The variational principle shows Z₂₄ ≥ 1. Suppose every |λ| were ≤ 0.97. Each λ⁶ would be ≤ 0.833, and the biggest Σλ²⁴ could be is 0.833⁴ + 0.233⁴ = 0.484 < 1. Contradiction.'],
            [26, 'So one eigenvalue λ₀ exceeds 0.97. Remove it: the others carry at most 1.0664 − 0.833 = 0.233 of the sixth moment, so Σ_rest |λ|²⁴ ≤ 0.233⁴ = 0.0030.'],
            [34, 'Step 3, the purity. Pₛ(24) ≥ (1 + 0.0030 / 0.97²⁴)⁻² = 0.988. A certified starting purity, defect ε ≈ 0.012, without ever knowing the spectrum (the true value is 0.99999).'],
            [41, 'The real certificate does exactly this with 12-site rings, 30-digit integers, symmetry sectors and a 72-site trial state, then runs six doublings to reach (2304, 784).']],
    draw(K, t) {
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
        [8.5, 'Z₆ = Σ λ⁶', '1.0664', 'exact, from small rings'],
        [10.5, 'max |λ|', '≤ 1.011', 'each λ⁶ ≤ Z₆'],
        [16.5, 'Z₂₄', '≥ 1', 'trial state (variational)'],
        [19, 'if all |λ| ≤ 0.97:  Σ λ²⁴', '≤ 0.484', 'fill the cap 0.97⁶ = 0.833 first'],
        [23, 'contradiction', '⇒ λ₀ > 0.97', ''],
        [27.5, 'rest: Σ λ⁶', '≤ 0.233', '1.0664 − 0.833'],
        [29.5, 'rest: Σ |λ|²⁴', '≤ 0.0030', '0.233⁴'],
        [35, 'Pₛ(24) ≥ (1 + 0.0030/0.97²⁴)⁻²', '= 0.988', 'ε ≈ 0.012'],
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
  });

  // 12 -----------------------------------------------------------------
  const DBLV = [[36, .0479158, .1815202], [72, .0604996, .1733812], [144, .0686606, .1632381], [288, .0713616, .1447092],
    [576, .0632939, .1055161], [1152, .0375793, .0445941], [2304, .0084513, .0054986]];
  chapter({
    id: 'seed', title: 'Where the starting numbers come from', dur: 22,
    beats: [[0, 'Small rings measure moments of the same eigenvalues: Zₖ = Σ λᵏ. Even powers are positive, so |λ| ≤ Zₖ^(1/k). Higher k squeezes the bound.'],
            [7, 'A trial state on 72 sites shows Z₇₂ > 1. No spectrum with every |λ| ≤ 0.999 can do that, so one eigenvalue dominates.'],
            [12, 'That gives purities at ring length 36. The doubling step itself, run six times, carries them to (2304, 784), both below 0.01.']],
    draw(K, t) {
      const {th} = K, A = 1 - seg(t, 11.5, 12.2);
      if (A > 0) {
        const lam = [1, .86, -.8, .62, .5, -.35, .2], k = 2 * Math.max(1, Math.round(lerp(1, 10, ease(seg(t, 1, 6.5))))), Zk = lam.reduce((s, l) => s + Math.pow(l, k), 0), bd = Math.min(1.7, Math.pow(Zk, 1 / k));
        const cx = 330, sc = 160, y = 170, X = v => cx + v * sc;
        K.line(X(-1.75), y, X(1.75), y, {color: th.mute, w: 1.5, alpha: A});
        [-1.5, -1, -.5, 0, .5, 1, 1.5].forEach(v => K.text(String(v), X(v), y + 32, {align: 'center', color: th.mute, size: 12, font: MONO, alpha: A}));
        K.rect(X(-bd), y - 18, X(bd) - X(-bd), 36, {fill: th.provedF, stroke: th.proved, w: 2, alpha: A});
        lam.forEach((l, i) => K.line(X(l), y - 14, X(l), y + 14, {color: i ? th.ink : th.space, w: i ? 2 : 3.5, alpha: A}));
        K.text(`k = ${k}:  |λ| ≤ ${Math.pow(Zk, 1 / k).toFixed(3)}`, X(1.8), y + 6, {font: MONO, size: 16, color: th.proved, alpha: A});
        K.text('hidden eigenvalues (ticks); the box is all a computer can certify from Zₖ', 60, 90, {color: th.mute, size: 14, alpha: A});
        const b = seg(t, 7, 7.8) * A;
        K.bars(80, 290, 300, 150, [.999, .999, .99, .98, .97, .96].map(v => v * .95), th.fail, {alpha: b, flat: true});
        K.bars(520, 290, 300, 150, [1, .6, .5, .45, .35, .3], th.ok, {alpha: b});
        K.text('every |λ| ≤ 0.999  ⇒  Z₇₂ ≤ 0.93', 80, 470, {color: th.fail, size: 15, alpha: b});
        K.text('trial state: Z₇₂ > 1  ⇒  one λ₀ > 0.999', 520, 470, {color: th.ok, size: 15, alpha: b});
      }
      const B = seg(t, 12.2, 13);
      if (B > 0) {
        const gx = 120, gy = 70, gw = 700, gh = 340, lo = -2.4, hi = -.6, X = j => gx + j / 6 * gw, Y = v => gy + gh * (hi - Math.log10(v)) / (hi - lo);
        K.axes(gx, gy, gw, gh, {alpha: B});
        [[.12731, 'ε* ≈ 0.127', th.mute], [.01, 'target 0.01', th.proved]].forEach(([v, nm, col]) => { K.line(gx, Y(v), gx + gw, Y(v), {color: col, dash: [6, 5], alpha: B, w: 1.5}); K.text(nm, gx + gw + 8, Y(v) + 5, {color: col, size: 14, alpha: B}); });
        DBLV.forEach((r, j) => K.text(String(r[0]), X(j), gy + gh + 22, {align: 'center', color: th.mute, size: 13, font: MONO, alpha: B}));
        K.text('ring length (β doubles alongside, 12¼ → 784)', gx, gy + gh + 46, {color: th.mute, size: 14, alpha: B});
        const prog = ease(seg(t, 13, 19)) * 6;
        [[1, th.space, 'spatial defect'], [2, th.time, 'thermal defect']].forEach(([k, col, nm], q) => {
          const pts = DBLV.map((r, j) => [X(j), Y(r[k])]).filter((_, j) => j <= prog);
          K.path(pts, {color: col, w: 3, alpha: B}); pts.forEach(p => K.dot(p[0], p[1], 5, col, B));
          K.text(nm, gx + 12, gy + 22 + 22 * q, {color: col, size: 15, alpha: B});
        });
      }
    }
  });

  // 13 -----------------------------------------------------------------
  chapter({
    id: 'open', title: 'Open chains', dur: 22,
    beats: [[0, 'An open spin-1 chain has a free spin-½ at each end. A small field on the end spins pins them.'],
            [6, 'A ring closes with a trace. An open chain closes with a boundary vector v. Apply the transfer matrix again and again …'],
            [11.5, '… the part of v along the top eigenvector ψ survives; the rest shrinks. Its squared length c is what the long chain sees.'],
            [16.5, 'If v missed ψ, nothing would survive. The open-chain proof tracks c and shows it stays positive at every scale.']],
    draw(K, t) {
      const {th} = K, A = 1 - seg(t, 5.6, 6.2);
      if (A > 0) {
        const N = 7, x0 = 220, dx = 86, y = 230, fon = seg(t, 3, 4);
        for (let i = 0; i < N; i++) {
          const x = x0 + i * dx;
          K.ctx.save(); K.ctx.globalAlpha = A; K.ctx.setLineDash([4, 3]); K.ctx.strokeStyle = th.mute; K.ctx.lineWidth = 1.5; K.ctx.beginPath(); K.ctx.ellipse(x, y, 30, 20, 0, 0, 7); K.ctx.stroke(); K.ctx.restore();
          if (i < N - 1) K.line(x + 13, y, x + dx - 13, y, {color: th.time, w: 3.5, alpha: A});
          K.dot(x - 13, y, 6.5, th.ink, A); K.dot(x + 13, y, 6.5, th.ink, A);
        }
        for (const x of [x0 - 13, x0 + (N - 1) * dx + 13]) { K.circle(x, y, 11, {color: th.proved, w: 3, alpha: A}); K.arrow(x, y - 18, x, y - 62, {color: th.space, w: 3, alpha: A * fon}); }
        K.text('each spin 1 = two spin-½; neighbours pair into singlets (blue)', W / 2, 320, {align: 'center', color: th.mute, size: 15, alpha: A});
        K.text('one spin-½ left at each end', W / 2, 350, {align: 'center', color: th.proved, size: 17, alpha: A * (1 - fon)});
        K.text('end fields pin them: one ground state', W / 2, 350, {align: 'center', color: th.space, size: 17, alpha: A * fon});
      }
      const B = seg(t, 6.2, 7);
      if (B > 0) {
        const ox = 330, oy = 380, sc = 260, ang = t < 16.5 ? 1.05 : Math.PI / 2;
        const sw = seg(t, 16.5, 17.2), aa = lerp(1.05, Math.PI / 2, sw);
        K.arrow(ox - 20, oy, ox + sc + 60, oy, {color: th.space, w: 2, alpha: B}); K.text('top eigenvector ψ', ox + sc + 64, oy + 5, {color: th.space, size: 15, alpha: B});
        K.arrow(ox, oy + 20, ox, oy - sc - 30, {color: th.mute, w: 2, alpha: B}); K.text('everything else', ox + 8, oy - sc - 34, {color: th.mute, size: 15, alpha: B});
        const vx = Math.cos(aa), vy = Math.sin(aa), steps = Math.floor(seg(t, 7.5, 11.5) * 6 + seg(t, 17.3, 20) * 6 * (t > 17.3));
        const r = .55; let px = vx, py = vy; const trail = [[px, py]];
        const nS = t < 16.5 ? Math.floor(seg(t, 7.5, 11.5) * 6) : Math.floor(seg(t, 17.3, 20.5) * 6);
        for (let k = 0; k < nS; k++) { py *= r; trail.push([px, py]); }
        trail.forEach(([a, b], k) => K.dot(ox + a * sc, oy - b * sc, k === trail.length - 1 ? 0 : 3.5, th.mute, B * .7));
        K.line(ox, oy, ox + vx * sc, oy - vy * sc, {color: th.proved, w: 1.5, dash: [5, 4], alpha: B * .5});
        const [ex, ey] = trail[trail.length - 1];
        K.arrow(ox, oy, ox + ex * sc, oy - ey * sc, {color: th.proved, w: 3.5, alpha: B});
        K.text(nS ? `after ${nS} more sites` : 'boundary vector v', ox + ex * sc + 12, oy - ey * sc - 8, {color: th.proved, size: 15, alpha: B});
        const c = seg(t, 12, 13) * (1 - sw);
        K.line(ox, oy + 2, ox + vx * sc, oy + 2, {color: th.space, w: 6, alpha: c * B});
        K.text(`√c = |⟨ψ, v⟩| = ${vx.toFixed(2)}`, ox + vx * sc / 2, oy + 34, {align: 'center', color: th.space, size: 16, font: MONO, alpha: c * B});
        K.text('the projection never shrinks', 640, 170, {color: th.ink, size: 17, alpha: seg(t, 13, 14) * (1 - sw)});
        K.text('v ⟂ ψ: the arrow shrinks to zero', 640, 170, {color: th.fail, size: 17, alpha: seg(t, 18, 19)});
      }
    }
  });

  // play the densest figures more slowly
  function slow(id, S) { const c = CH.find(x => x.id === id); if (!c) return; const d0 = c.draw; c.dur *= S; c.beats = c.beats.map(([b, txt]) => [b * S, txt]); c.draw = (K, t) => d0(K, t / S); }
  slow('rptest', 1.35); slow('seed', 1.35);
  root.STBV = {CH, kit, W, H, minusH, realign, eigSym, slow};
})(typeof window !== 'undefined' ? window : globalThis);
