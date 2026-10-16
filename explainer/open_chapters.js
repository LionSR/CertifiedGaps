// Part II: open chains with boundary fields (Appendices A–C of notes/thermal-mpo-gap-note.tex).
// Same drawing kit and 960x540 canvas as chapters.js; these figures form their own film.
(function (root) {
  'use strict';
  const {W, seg, ease, lerp, slow} = root.STBV;
  const OPEN = [];
  const chapter = c => OPEN.push(c);

  // valence-bond cartoon of an open spin-1 chain
  function vbChain(K, x0, y, n, dx, alpha, fieldOn) {
    const {th, ctx} = K;
    for (let i = 0; i < n; i++) {
      const x = x0 + i * dx;
      ctx.save(); ctx.globalAlpha *= alpha; ctx.setLineDash([3, 3]); ctx.strokeStyle = th.mute; ctx.lineWidth = 1; ctx.beginPath(); ctx.ellipse(x, y, dx * .34, 17, 0, 0, 7); ctx.stroke(); ctx.restore();
      if (i < n - 1) K.line(x + 12, y, x + dx - 12, y, {color: th.time, w: 3, alpha});
      K.dot(x - 12, y, 6, th.ink, alpha); K.dot(x + 12, y, 6, th.ink, alpha);
    }
    for (const x of [x0 - 12, x0 + (n - 1) * dx + 12]) {
      K.circle(x, y, 10, {color: th.proved, w: 2.5, alpha});
      if (fieldOn > 0) K.arrow(x, y - 16, x, y - 56, {color: th.space, w: 3, alpha: alpha * fieldOn});
    }
  }
  // a strip of transfer boxes between two caps; returns box centers
  function strip(K, x0, y, n, dx, o = {}) {
    const {th} = K, al = o.alpha ?? 1, xs = [];
    K.line(x0 - dx * .6, y, x0 + (n - 1) * dx + dx * .6, y, {color: th.mute, w: 1.5, alpha: al});
    for (let i = 0; i < n; i++) {
      const x = x0 + i * dx; xs.push(x);
      const hl = o.hl ? o.hl(i) : null;
      K.rect(x - 15, y - 15, 30, 30, {fill: hl || th.bg, stroke: o.color || th.space, w: 2, alpha: al});
      K.text(o.label ? o.label(i) : 'X', x, y + 6, {align: 'center', color: o.color || th.space, size: 15, alpha: al});
    }
    if (o.caps !== false) {
      for (const [x, l] of [[x0 - dx, o.left || 'v̄'], [x0 + n * dx, o.right || 'v']]) {
        K.circle(x, y, 15, {color: th.proved, w: 2.5, alpha: al});
        K.text(l, x, y + 6, {align: 'center', color: th.proved, size: 16, alpha: al});
      }
    }
    return xs;
  }

  // 1 -------------------------------------------------------------------
  const EDLV = {"0.0": {"4": [0.0, 0.5092, 0.5092, 0.5092], "5": [0.0, 0.0, 0.0, 0.5466], "6": [0.0, 0.3078, 0.3078, 0.3078], "7": [0.0, 0.0, 0.0, 0.331], "8": [0.0, 0.2019, 0.2019, 0.2019], "9": [0.0, 0.0, 0.0, 0.2127], "10": [0.0, 0.1383, 0.1383, 0.1383], "11": [0.0, 0.0, 0.0, 0.1418], "12": [0.0, 0.0971, 0.0971, 0.0971]}, "0.6": {"4": [0.0, 0.0707, 0.6484, 1.0177], "5": [0.0, 0.7187, 1.2976, 1.3842], "6": [0.0, 0.2496, 0.6512, 0.8365], "7": [0.0, 0.6365, 1.0479, 1.081], "8": [0.0, 0.3354, 0.6425, 0.7332], "9": [0.0, 0.5832, 0.8965, 0.8988], "10": [0.0, 0.3815, 0.6296, 0.6668], "11": [0.0, 0.5461, 0.7802, 0.7967], "12": [0.0, 0.4075, 0.6154, 0.6205]}}; // open spin-1 chain, lowest levels above the ground state (scipy ED, floating point)
  chapter({
    id: 'oedge', title: 'Ends carry free spins', dur: 26,
    beats: [[0, 'Cut the ring open. In the valence-bond picture each spin 1 is two spin-½, and neighbours pair into singlets (blue). One spin-½ is left over at each end.'],
            [5.5, 'Real numbers, no field (exact diagonalization, n = 4 … 12). Even n: three states sit just above the ground state, at 0.51, 0.31, 0.20, 0.14, 0.10. Odd n: the ground state itself is three-fold degenerate. No uniform gap.'],
            [13.5, 'Now put h = 3/5 on the two end spins. For odd n the ground state is alone, with gap 0.72, 0.64, 0.58, 0.55, levelling off. Even n are worse: the two ends prefer opposite directions.'],
            [20.5, 'That is the chain of the theorem: H = Σ S·S − h(S₁ᶻ + Sₙᶻ), h = 3/5, n odd.']],
    draw(K, t) {
      const {th} = K, fon = seg(t, 13.5, 14.5);
      vbChain(K, 180, 90, 7, 100, seg(t, 0, .8), fon);
      const gx = 160, gy = 180, gw = 560, gh = 260, X = n => gx + (n - 4) / 8 * gw, Y = e => gy + gh - e / 1.1 * gh, b = seg(t, 5.5, 6.3);
      K.line(gx, gy + gh, gx + gw, gy + gh, {color: th.mute, w: 1.5, alpha: b}); K.line(gx, gy, gx, gy + gh, {color: th.mute, w: 1.5, alpha: b});
      for (let n = 4; n <= 12; n++) K.text(String(n), X(n), gy + gh + 22, {align: 'center', color: n % 2 ? th.ink : th.mute, size: 14, alpha: b});
      [0, .5, 1].forEach(v => K.text(v.toFixed(1), gx - 8, Y(v) + 5, {align: 'right', color: th.mute, size: 13, alpha: b}));
      K.text('energy above the ground state', gx, gy - 10, {color: th.mute, size: 13, alpha: b});
      K.text('chain length n', gx + gw, gy + gh + 44, {align: 'right', color: th.mute, size: 13, alpha: b});
      const plot = (key, col, al, shown) => {
        for (let n = 4; n <= 12; n++) {
          if (n - 4 > shown) continue;
          const L = EDLV[key][n];
          for (let k = 1; k < 4; k++) { const e = L[k]; if (e > 1.05) continue; K.line(X(n) - 16, Y(e), X(n) + 16, Y(e), {color: col, w: 3, alpha: al}); }
          K.line(X(n) - 16, Y(0), X(n) + 16, Y(0), {color: th.ink, w: 3, alpha: al});
        }
      };
      const s0 = (t - 6) * 2.2, s1 = (t - 14) * 2.2;
      plot('0.0', th.time, b * (1 - .75 * fon), s0);
      if (fon > 0) {
        plot('0.6', th.space, fon, s1);
        const odd = [5, 7, 9, 11];
        odd.forEach((n, i) => { if (n - 4 <= s1) K.text(EDLV['0.6'][n][1].toFixed(2), X(n), Y(EDLV['0.6'][n][1]) - 10, {align: 'center', color: th.space, size: 13, alpha: fon}); });
      }
      K.text(fon < .5 ? 'no field' : 'end fields h = 3/5', gx + gw + 14, gy + 20, {color: fon < .5 ? th.time : th.space, size: 16, alpha: b});
      K.text('lines: the three lowest excited levels', gx + gw + 14, gy + 44, {color: th.mute, size: 13, alpha: b});
      K.text('floating-point exact diagonalization', gx + gw + 14, gy + 64, {color: th.mute, size: 13, italic: true, alpha: b});
      K.tex('H_n=\\sum_{i=1}^{n-1}\\mathbf S_i\\cdot\\mathbf S_{i+1}-h\\,(S_1^z+S_n^z),\\qquad h=\\tfrac35', W / 2, 505, {align: 'center', size: 22, alpha: seg(t, 20.5, 21.5)});
    }
  });

  // 2 -------------------------------------------------------------------
  chapter({
    id: 'ocyl', title: 'A cylinder, not a torus', dur: 17,
    beats: [[0, 'For a ring, space and imaginary time are both circles: the torus.'],
            [4, 'An open chain is still periodic in imaginary time, because Z is still a trace. But space now has two ends: the spacetime is a cylinder.'],
            [10, 'Each end of the cylinder is closed by a boundary vector. A ring closes with a trace; an open chain closes with ⟨v̄| X X ⋯ X |v⟩.']],
    draw(K, t) {
      const {th, ctx} = K, m = 2 - ease(seg(t, 4, 8));
      STB.drawTorus(ctx, {cx: W / 2, cy: 220, scale: lerp(240, 400, seg(m, 1, 2)), m, yaw: -0.5 + .1 * Math.sin(t * .5), pitch: lerp(.55, .9, seg(m, 1, 2)),
        cols: 16, rows: 8, colHi: -1, nu: 64, nv: 32, lineW: 1.2, shadeAmt: .25, colors: {space: th.space, time: th.time, fill: th.torusFill, hi: th.proved}});
      K.text(m > 1.5 ? 'ring: torus' : 'open chain: cylinder', W / 2, 410, {align: 'center', size: 18, alpha: 1});
      const b = seg(t, 10, 11);
      strip(K, 290, 470, 7, 64, {alpha: b});
      K.text('boundary', 290 - 64, 506, {align: 'center', color: th.proved, size: 13, alpha: b});
      K.text('boundary', 290 + 7 * 64, 506, {align: 'center', color: th.proved, size: 13, alpha: b});
    }
  });

  // 3 -------------------------------------------------------------------
  chapter({
    id: 'ovec', title: 'What a boundary vector is', dur: 23,
    beats: [[0, 'The end field acts on one site only. Expanding e^(−βH) puts field events on the end site’s world line; each inserts G_z = iSᶻ.'],
            [6, 'The list of field events is a history a, like a bond’s. Its weight is the boundary vector: v(a) = (ih)ᵈ when all d events are of type z.'],
            [12, 'With field −h at one end and +h at the other, the two caps are the same vector, one conjugated. The chain is ⟨v̄|Xⁿ|v⟩, and its even powers are squared norms.'],
            [18, 'The equal-field chain of the theorem is reached by a π rotation U on one cap, which flips the sign of that end field.']],
    draw(K, t) {
      const {th} = K, A = 1 - seg(t, 11.6, 12.2);
      if (A > 0) {
        const hist = [[], [.5], [.25, .7], [.2, .45, .8]], y0 = 90, y1 = 380;
        hist.forEach((ev, k) => {
          const x = 140 + k * 200, al = A * seg(t, .5 + k * .6, 1.1 + k * .6);
          K.line(x - 46, y0, x - 46, y1, {color: th.ink, w: 3, alpha: al});
          K.line(x - 52, y0, x - 52, y1, {color: th.ink, w: 1, alpha: al});
          K.line(x, y0, x, y1, {color: th.mute, w: 2, alpha: al});
          ev.forEach(f => { const y = y1 - f * (y1 - y0); K.line(x - 46, y, x, y, {color: th.space, w: 2.5, alpha: al}); K.dot(x, y, 5, th.space, al); K.text('G_z', x + 8, y + 5, {color: th.ink, size: 13, alpha: al}); });
          K.text(`d = ${ev.length}`, x - 23, y1 + 24, {align: 'center', color: th.mute, size: 14, alpha: al});
          const val = ['1', 'ih', '(ih)² = −h²', '(ih)³ = −ih³'][k];
          K.text(`v(a) = ${val}`, x - 23, y1 + 50, {align: 'center', color: th.proved, size: 15, alpha: al * seg(t, 6.5, 7.5)});
        });
        K.text('field', 140 - 56, y0 - 12, {align: 'right', color: th.ink, size: 13, alpha: A});
        K.text('end site', 140 + 4, y0 - 12, {color: th.mute, size: 13, alpha: A});
        K.text('τ ↑', 50, y0 + 10, {color: th.time, size: 14, alpha: A});
      }
      const B = seg(t, 12.2, 12.9);
      if (B > 0) {
        strip(K, 300, 200, 6, 70, {alpha: B});
        K.text('field −h', 300 - 70, 150, {align: 'center', color: th.space, size: 15, alpha: B});
        const flip = seg(t, 18.5, 20);
        K.text(flip < .5 ? 'field +h' : 'field −h', 300 + 6 * 70, 150, {align: 'center', color: th.space, size: 15, alpha: B});
        K.arrow(300 - 70, 172, 300 - 70, 132 + 6, {color: th.space, w: 2.5, alpha: B * 0});
        K.text('same vector, conjugated: every even power is ⟨Φ|Φ⟩ ≥ 0', W / 2, 270, {align: 'center', color: th.ok, size: 16, alpha: B * seg(t, 14, 15)});
        const ux = 300 + 6 * 70 + 50;
        K.rect(ux - 16, 184, 32, 32, {stroke: th.time, w: 2, alpha: B * seg(t, 18, 18.8)});
        K.text('U', ux, 206, {align: 'center', color: th.time, size: 14, alpha: B * seg(t, 18, 18.8)});
        K.text('π rotation: flips the right field, so both ends have −h', W / 2, 330, {align: 'center', color: th.time, size: 16, alpha: B * seg(t, 19, 20)});
        K.text('opposite fields:  ⟨v̄| Xⁿ |v⟩        equal fields:  ⟨v̄| Xⁿ U |v⟩', W / 2, 420, {align: 'center', size: 16, alpha: B * seg(t, 20, 21)});
      }
    }
  });

  // 4 -------------------------------------------------------------------
  chapter({
    id: 'over', title: 'The overlap', dur: 23,
    beats: [[0, 'Split the normalized transfer matrix: T = |ψ⟩⟨ψ| + R, with ψ the top eigenvector and ‖R‖ < 1.'],
            [5, 'Apply T to the boundary vector, again and again. The part along ψ stays; the rest shrinks by ‖R‖ each time.'],
            [11, 'So ⟨v|Tᵐ|v⟩ = c + Dₘ: a fixed overlap c = |⟨ψ, v⟩|² plus a remainder that dies away with the length m.'],
            [17, 'The boundary defect e = Dₘ / Bₘ is the share not carried by ψ. Any e < 1 proves c > 0: a long chain does see the ground state.']],
    draw(K, t) {
      const {th} = K, ox = 120, oy = 400, sc = 280, vx = Math.cos(1.05), vy = Math.sin(1.05);
      K.arrow(ox - 20, oy, ox + sc + 40, oy, {color: th.space, w: 2}); K.text('ψ', ox + sc + 48, oy + 6, {color: th.space, size: 20});
      K.arrow(ox, oy + 20, ox, oy - sc - 10, {color: th.mute, w: 2}); K.text('everything else', ox + 8, oy - sc - 14, {color: th.mute, size: 14});
      const nS = Math.floor(seg(t, 5.5, 10.5) * 7), trail = [[vx, vy]];
      for (let k = 0; k < nS; k++) trail.push([vx, trail[trail.length - 1][1] * .6]);
      trail.forEach(([a, b], k) => K.dot(ox + a * sc, oy - b * sc, 3.5, th.mute, k < trail.length - 1 ? .8 : 0));
      K.line(ox, oy, ox + vx * sc, oy - vy * sc, {color: th.proved, dash: [5, 4], w: 1.5, alpha: nS ? .5 : 0});
      const [ex, ey] = trail[trail.length - 1];
      K.arrow(ox, oy, ox + ex * sc, oy - ey * sc, {color: th.proved, w: 3.5});
      K.text(nS ? `Tᵐv, m = ${nS}` : 'v', ox + ex * sc + 10, oy - ey * sc - 8, {color: th.proved, size: 16});
      K.line(ox, oy + 3, ox + vx * sc, oy + 3, {color: th.space, w: 7, alpha: seg(t, 11, 12)});
      K.text('√c', ox + vx * sc / 2, oy + 26, {align: 'center', color: th.space, size: 17, alpha: seg(t, 11, 12)});
      const b = seg(t, 11, 11.8);
      if (b > 0) {
        const gx = 520, gy = 110, gw = 380, gh = 270, Y = v => gy + gh - v / 1.05 * gh, c = vx * vx;
        K.line(gx, gy + gh, gx + gw, gy + gh, {color: th.mute, w: 1.5, alpha: b});
        for (let m = 0; m <= 14; m += 2) {
          const D = vy * vy * Math.pow(.6, m), x = gx + 10 + m / 14 * (gw - 30);
          K.rect(x, Y(c), 14, Y(0) - Y(c), {fill: th.space, alpha: b * .75});
          K.rect(x, Y(c + D), 14, Y(c) - Y(c + D), {fill: th.mute, alpha: b * .5});
          K.text(String(m), x + 7, gy + gh + 18, {align: 'center', color: th.mute, size: 12, alpha: b});
        }
        K.text('chain length m', gx + gw, gy + gh + 40, {align: 'right', color: th.mute, size: 13, alpha: b});
        K.text('overlap c', gx + gw - 4, Y(c / 2) + 5, {align: 'right', color: th.bg, size: 14, alpha: 0});
        K.text('c (red) stays', gx, gy - 28, {color: th.space, size: 15, alpha: b});
        K.text('remainder Dₘ (grey) dies', gx, gy - 8, {color: th.mute, size: 15, alpha: b});
        K.text('e = Dₘ / (c + Dₘ)  <  1   ⇒   c > 0', gx, gy + gh + 72, {color: th.proved, size: 17, alpha: seg(t, 17, 18)});
        const tb = seg(t, 18.5, 19.3);
        if (tb > 0) { const rows = [0, 2, 4, 8, 16];
          rows.forEach((m, i) => { const D = .75 * Math.pow(.6, m), B = .25 + D, yy = 494 + 0;
            const xx = 60 + i * 175; K.text(`m = ${m}:  e = ${(D / B).toFixed(3)}`, xx, 500, {size: 15, alpha: tb * seg(t, 19 + i * .5, 19.5 + i * .5)}); });
          K.text('worked example: c = 0.25, remainder 0.75 · 0.6ᵐ', 60, 528, {color: th.mute, size: 13, italic: true, alpha: tb}); }
      }
    }
  });

  // 5 -------------------------------------------------------------------
  chapter({
    id: 'ofold', title: 'Why the remainder is a square', dur: 21,
    beats: [[0, 'The remainder is Dₘ = ⟨v|Rᵐ|v⟩ for even m. Fold the strip at its middle.'],
            [5.5, 'The two halves are mirror images, the reflection-positivity argument once more: Dₘ = ‖R^(m/2) v‖² ≥ 0.'],
            [11, 'Now absorb k/2 factors into each end, w = R^(k/2) v. The middle has norm at most ‖R‖^(m−k), so Dₘ ≤ ‖R‖^(m−k) Dₖ.'],
            [16.5, 'Longer chains shrink the remainder geometrically, while the overlap c stays put. That is what lets the boundary defect improve with length.']],
    draw(K, t) {
      const {th} = K, n = 8, dx = 66, x0 = W / 2 - (n - 1) * dx / 2, y = 250, mid = x0 + (n - 1) * dx / 2;
      const f = ease(seg(t, 2, 5.5)), C = 1 - seg(t, 10.6, 11.2);
      if (C > 0) {
        for (let i = 0; i < n; i++) {
          const x = x0 + i * dx, right = i >= n / 2;
          const xx = right ? lerp(x, 2 * mid - x, f) : x, yy = right ? lerp(y, y - 90, f) : y;
          K.rect(xx - 15, yy - 15, 30, 30, {fill: th.bg, stroke: th.space, w: 2, alpha: C});
          K.text('R', xx, yy + 6, {align: 'center', color: th.space, size: 15, alpha: C});
        }
        const lx = x0 - dx, rx = x0 + n * dx, rxx = lerp(rx, 2 * mid - rx, f), ryy = lerp(y, y - 90, f);
        K.circle(lx, y, 15, {color: th.proved, w: 2.5, alpha: C}); K.text('v̄', lx, y + 6, {align: 'center', color: th.proved, size: 16, alpha: C});
        K.circle(rxx, ryy, 15, {color: th.proved, w: 2.5, alpha: C}); K.text('v', rxx, ryy + 6, {align: 'center', color: th.proved, size: 16, alpha: C});
        K.line(mid, y - 120, mid, y + 40, {color: th.mute, dash: [5, 4], w: 1.5, alpha: C * (1 - f)});
        if (f > .98) { K.line(x0 + (n / 2 - 1) * dx + 20, y, x0 + (n / 2 - 1) * dx + 46, y - 45, {color: th.mute, w: 1.5, alpha: C}); K.line(x0 + (n / 2 - 1) * dx + 46, y - 45, x0 + (n / 2 - 1) * dx + 20, y - 90, {color: th.mute, w: 1.5, alpha: C}); }
        K.text('fold', mid, y + 60, {align: 'center', color: th.mute, size: 14, alpha: C * (1 - f)});
        K.text('top row = bottom row mirrored:   Dₘ = ‖R^(m/2) v‖² ≥ 0', W / 2, 400, {align: 'center', color: th.ok, size: 18, alpha: C * seg(t, 6, 7)});
      }
      const D = seg(t, 11.2, 12);
      if (D > 0) {
        const m = 10, k = 4, d2 = 62, xs0 = W / 2 - (m - 1) * d2 / 2;
        strip(K, xs0, 220, m, d2, {alpha: D, label: () => 'R', hl: i => (i < k / 2 || i >= m - k / 2) ? th.provedF : null});
        K.rect(xs0 - d2 - 22, 180, d2 * (k / 2 + 1) + 14, 80, {stroke: th.proved, w: 1.5, dash: [5, 4], alpha: D * seg(t, 12, 13)});
        K.rect(xs0 + (m - k / 2) * d2 - 22, 180, d2 * (k / 2 + 1) + 14, 80, {stroke: th.proved, w: 1.5, dash: [5, 4], alpha: D * seg(t, 12, 13)});
        K.text('w', xs0 - 10, 290, {align: 'center', color: th.proved, size: 20, alpha: D * seg(t, 12, 13)});
        K.text('w', xs0 + (m - 1) * d2 + 10, 290, {align: 'center', color: th.proved, size: 20, alpha: D * seg(t, 12, 13)});
        K.text('middle: norm ≤ ‖R‖^(m−k)', W / 2, 160, {align: 'center', color: th.space, size: 16, alpha: D * seg(t, 13, 14)});
        K.text('Dₘ ≤ ‖R‖^(m−k) · Dₖ', W / 2, 360, {align: 'center', size: 24, alpha: D * seg(t, 14, 15)});
        const gx = 330, gy = 400, gw = 300, gh = 80;
        for (let mm = 4; mm <= 20; mm += 2) { const v = Math.pow(.7, mm - 4), x = gx + (mm - 4) / 16 * gw; K.rect(x, gy + gh - v * gh, 10, v * gh, {fill: th.mute, alpha: D * seg(t, 16.5, 17.5) * .6}); }
        K.text('remainder versus m', gx + gw + 14, gy + gh, {color: th.mute, size: 13, alpha: D * seg(t, 16.5, 17.5)});
      }
    }
  });

  // 6 -------------------------------------------------------------------
  chapter({
    id: 'omoves', title: 'Three moves per scale', dur: 18,
    beats: [[0, 'At each scale j track two numbers: the boundary defect e, and the open-chain thermal purity at length Mⱼ. Both start at most x.'],
            [3.5, 'Move 1: the overlap recurrence. The longer chain at the next temperature squeezes the boundary defect to at most 10x².'],
            [7.5, 'Move 2: double the length. With the new overlap bound, the open-chain purity at 2Mⱼ costs only about 2x.'],
            [11, 'Move 3: square the Gibbs state, which doubles β. The defect drops to 10x².'],
            [14, 'Out comes the same pair at the next scale with x replaced by 10x². The step repeats, as on the ring.']],
    draw(K, t) {
      const {th} = K, xe = 90, xp = 520, bw = 350, bh = 50, rows = [80, 175, 270, 365];
      const box = (x, r, col, s, at, dash) => { const al = seg(t, at, at + .7); K.rect(x, rows[r], bw, bh, {fill: th.bg, stroke: col, w: 2, dash, alpha: al}); K.text(s, x + 16, rows[r] + 32, {size: 16, alpha: al}); };
      box(xe, 0, th.mute, 'boundary defect, scale j:  e ≤ x', .3); box(xp, 0, th.time, 'open-chain purity at (Mⱼ, βⱼ):  ≥ 1 − x', .3);
      box(xe, 1, th.mute, 'boundary defect, scale j+1:  e ≤ 10x²', 4);
      box(xp, 2, th.time, 'purity at (2Mⱼ, βⱼ):  ≥ 1 − 2x − 16x²', 8, [6, 5]);
      box(xp, 3, th.time, 'purity at (2Mⱼ, βⱼ₊₁):  ≥ 1 − 10x²', 11.5);
      const mv = (x1, y1, x2, y2, k, lab, at) => { const al = seg(t, at - .8, at); K.arrow(x1, y1, x2, y2, {color: th.proved, w: 2.5, alpha: al}); K.text(`${k}  ${lab}`, (x1 + x2) / 2 + 12, (y1 + y2) / 2 + 5, {color: th.proved, size: 15, alpha: al}); };
      mv(xe + 40, rows[0] + bh, xe + 40, rows[1] - 2, '①', 'overlap recurrence', 4);
      mv(xp + 40, rows[0] + bh, xp + 40, rows[2] - 2, '②', 'double the length', 8);
      mv(xp + 40, rows[2] + bh, xp + 40, rows[3] - 2, '③', 'square in time', 11.5);
      K.arrow(xp, rows[0] + bh / 2, xe + bw - 20, rows[1] + 6, {color: th.mute, w: 1.5, alpha: seg(t, 3.5, 4)});
      K.arrow(xe + bw, rows[1] + bh / 2, xp - 6, rows[2] + 12, {color: th.mute, w: 1.5, alpha: seg(t, 7.5, 8)});
      K.text('(M, β, x)  →  (2M, 2β, 10x²)', xe, 470, {size: 24, color: th.proved, alpha: seg(t, 14, 15)});
    }
  });

  // 6b ------------------------------------------------------------------
  chapter({
    id: 'onum', title: 'One open-chain step, in numbers', dur: 30,
    beats: [[0, 'The note’s first step, j = 2 → 3, with every number. Given: boundary defect E = 0.14, open-chain thermal defect η = 0.21, periodic defects ε₂ = 0.052 and ε₃ = 6ε₂² = 0.0162.'],
            [8, 'Move 1. The suppression factor is q = 0.142. The overlap recurrence gives E′ ≤ q/(1−q) · [1/((1−E)²(1−η)) − 1] = 0.1655 × 0.7115 = 0.118, below the table entry 0.12.'],
            [17, 'Moves 2 and 3. With E′ rounded up to 0.12, the purity after doubling is at least (1−η)(1−ε₂)(1−ε₃)(1−E′) = 0.79 × 0.948 × 0.984 × 0.88 = 0.648; squaring turns the defect 0.352 into f(0.352) = 0.147, below 0.155.'],
            [25, 'Both defects went down. Two more such steps reach 0.001; after that the clean rule x → 10x² takes over.']],
    draw(K, t) {
      const {th} = K;
      const lines = [
        [1, 'given', 'E = 0.14     η = 0.21     ε₂ = 0.052     ε₃ = 0.0162', th.ink],
        [8.5, 'move 1', 'q = 0.142', th.ink],
        [10.5, '', 'E′ ≤ q/(1−q) · [ 1 / ((1−E)²(1−η)) − 1 ]', th.ink],
        [12.5, '', '   = 0.1655 × [ 1 / (0.7396 × 0.79) − 1 ] = 0.1655 × 0.7115 = 0.118', th.proved],
        [17.5, 'move 2', '(1−η)(1−ε₂)(1−ε₃)(1−0.12) = 0.79 × 0.948 × 0.984 × 0.88 = 0.648', th.ink],
        [20.5, 'move 3', 'η′ ≤ f(1 − 0.648) = f(0.352) = 0.352² / (2 · 0.648²) = 0.147', th.proved],
        [25.5, 'result', 'E: 0.14 → 0.118 ≤ 0.12        η: 0.21 → 0.147 ≤ 0.155', th.ok],
      ];
      lines.forEach(([at, tag, txt, col], i) => { const al = seg(t, at, at + .8), y = 110 + i * 52;
        K.text(tag, 80, y, {color: th.mute, size: 15, italic: true, alpha: al}); K.text(txt, 180, y, {color: col, size: 18, alpha: al}); });
    }
  });

  // 7 -------------------------------------------------------------------
  chapter({
    id: 'ostart', title: 'Starting numbers', dur: 19,
    beats: [[0, 'The open-chain iteration starts at β = 49 from four certified numbers: two periodic purities at length 120 and two open-chain bounds from the boundary-field paper.'],
            [6.5, 'Three rounds of the moves take the boundary defect from 0.14 and the thermal defect from 0.21 down to 0.001 at β = 392, length 1920.'],
            [12.5, 'From there x → 10x²: 10⁻², 10⁻³, 10⁻⁵, 10⁻⁹, 10⁻¹⁷, doubly exponential, exactly as on the ring.']],
    draw(K, t) {
      const {th} = K, A = seg(t, 0, .8);
      const rows = [['Pₛ(120, 49) ≥ 0.948', 'periodic'], ['Pτ(240, 49) ≥ 0.948', 'periodic'], ['e₁₄₂ ≤ 0.14', 'open chain'], ['open-chain Pτ(240, 49) ≥ 0.79', 'open chain']];
      rows.forEach(([s, k], i) => { const al = A * seg(t, .5 + i * .6, 1.1 + i * .6); K.text(s, 70, 110 + i * 40, {size: 17, alpha: al}); K.text(k, 330, 110 + i * 40, {color: th.mute, size: 14, italic: true, alpha: al}); });
      const E = [[2, .14, .21], [3, .12, .155], [4, .016, .026], [5, .001, .001]], XI = [];
      for (let j = 5; j <= 9; j++) XI.push([j, Math.pow(10, -1 - Math.pow(2, j - 5))]);
      const gx = 470, gy = 80, gw = 420, gh = 330, lo = -18, X = j => gx + (j - 2) / 7 * gw, Y = v => gy + Math.log10(v) / lo * gh, b = seg(t, 6.5, 7.3);
      K.line(gx, gy + gh, gx + gw, gy + gh, {color: th.mute, w: 1.5, alpha: b}); K.line(gx, gy, gx, gy + gh, {color: th.mute, w: 1.5, alpha: b});
      [0, -5, -10, -15].forEach(v => K.text(v ? `10${String(v).split('').map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(c)]).join('')}` : '1', gx - 8, Y(Math.pow(10, v)) + 5, {align: 'right', color: th.mute, size: 13, alpha: b}));
      const upto = 2 + (seg(t, 7, 12) * 3) + seg(t, 13, 18) * 4.01;
      [[1, th.mute, 'table: boundary defect E'], [2, th.time, 'table: thermal defect η']].forEach(([k, col, nm], q) => {
        const pts = E.filter(r => r[0] <= upto).map(r => [X(r[0]), Y(r[k])]);
        K.path(pts, {color: col, w: 2.5, alpha: b, dash: q ? null : [6, 4]}); pts.forEach(p => K.dot(p[0], p[1], 4.5, col, b));
        K.text(nm, gx + gw - 4, gy + 20 + 20 * q, {align: 'right', color: col, size: 14, alpha: b});
      });
      const xp = XI.filter(r => r[0] <= Math.max(5, upto - 0)).map(r => [X(r[0]), Y(r[1])]), c2 = b * seg(t, 12.5, 13.3);
      K.path(xp, {color: th.ink, w: 2.5, alpha: c2}); xp.forEach(p => K.dot(p[0], p[1], 4.5, th.ink, c2));
      K.text('induction bound xⱼ, x₅ = 0.01', X(5.15), Y(1e-14), {color: th.ink, size: 14, alpha: c2});
      [2, 3, 4, 5, 6, 7, 8, 9].forEach(j => K.text(`j=${j}`, X(j), gy + gh + 20, {align: 'center', color: th.mute, size: 12, alpha: b}));
      K.line(X(5), gy, X(5), gy + gh, {color: th.proved, dash: [4, 4], w: 1.5, alpha: b * seg(t, 11, 12)});
      K.text('β = 392, M = 1920', X(5) + 6, gy + gh - 10, {color: th.proved, size: 13, alpha: b * seg(t, 11, 12)});
      K.text('table of the note (j ≤ 5); the induction restarts from x₅ = 0.01, above the table values', gx, gy + gh + 44, {color: th.mute, size: 13, alpha: b});
    }
  });

  // 8 -------------------------------------------------------------------
  chapter({
    id: 'ogap', title: 'Odd lengths, and the gap', dur: 24,
    beats: [[0, 'The theorem is for odd n, where Rⁿ need not be positive. Put M/2 factors of R on each end: w = R^(M/2) v.'],
            [6, 'What is left in the middle, R^(n−M) U, has norm at most 1. So the remainder is at most ‖w‖² = D(M), whatever its sign.'],
            [11.5, 'Every n ≥ 1920 lies in a window [Mⱼ, 2Mⱼ] with βⱼ = 392·2^(j−5). In that window the thermal purity is at least 1 − 4xⱼ > ½.'],
            [18, 'So the ground state is unique and e^(−βⱼΔ) < 10^(−2^(j−5)). Dividing by βⱼ: Δₙ > log 10 / 392 ≈ 0.0059 for every odd n ≥ 1921.']],
    draw(K, t) {
      const {th} = K, A = 1 - seg(t, 11.2, 11.8);
      if (A > 0) {
        const n = 9, k = 4, d2 = 62, xs0 = W / 2 - (n - 1) * d2 / 2;
        strip(K, xs0, 220, n, d2, {alpha: A, label: i => i === n - k / 2 - 1 + 1 - 1 && false ? 'U' : 'R', hl: i => (i < k / 2 || i >= n - k / 2) ? th.provedF : null});
        K.text(`n = ${n} (odd)`, W / 2, 140, {align: 'center', color: th.mute, size: 15, alpha: A});
        const g = seg(t, 1.5, 2.5);
        K.rect(xs0 - d2 - 22, 180, d2 * (k / 2 + 1) + 14, 80, {stroke: th.proved, w: 1.5, dash: [5, 4], alpha: A * g});
        K.rect(xs0 + (n - k / 2) * d2 - 22, 180, d2 * (k / 2 + 1) + 14, 80, {stroke: th.proved, w: 1.5, dash: [5, 4], alpha: A * g});
        K.text('w', xs0 - 10, 290, {align: 'center', color: th.proved, size: 20, alpha: A * g});
        K.text('w', xs0 + (n - 1) * d2 + 10, 290, {align: 'center', color: th.proved, size: 20, alpha: A * g});
        const ux = xs0 + (n - k / 2 - .5) * d2, u = seg(t, 6, 7);
        K.rect(ux - 11, 209, 22, 22, {fill: th.bg, stroke: th.time, w: 1.5, alpha: A * u}); K.text('U', ux, 226, {align: 'center', color: th.time, size: 13, alpha: A * u});
        K.text('middle: ‖R^(n−M) U‖ ≤ 1', W / 2, 330, {align: 'center', color: th.space, size: 17, alpha: A * u});
        K.text('|remainder| ≤ ‖w‖² = D(M)', W / 2, 380, {align: 'center', size: 24, alpha: A * seg(t, 7.5, 8.5)});
      }
      const B = seg(t, 11.8, 12.6);
      if (B > 0) {
        const lx = 80, lw = 800, ly = 230, lo = Math.log10(1500), hi = Math.log10(50000), X = v => lx + (Math.log10(v) - lo) / (hi - lo) * lw;
        K.line(lx, ly, lx + lw, ly, {color: th.mute, w: 1.5, alpha: B});
        [1920, 3840, 7680, 15360, 30720].forEach((M, i) => {
          const al = B * seg(t, 12.2 + i * .7, 12.8 + i * .7), j = 5 + i;
          K.rect(X(M), ly - 28, X(Math.min(2 * M, 50000)) - X(M), 28, {fill: i % 2 ? th.timeF : th.spaceF, alpha: al});
          K.line(X(M), ly - 34, X(M), ly + 6, {color: th.ink, w: 1.5, alpha: al});
          K.text(String(M), X(M), ly + 24, {align: 'center', color: th.mute, size: 13, alpha: al});
          K.text(`j = ${j}`, (X(M) + X(Math.min(2 * M, 50000))) / 2, ly - 40, {align: 'center', size: 14, alpha: al});
          K.text(`β = ${392 * Math.pow(2, i)}`, (X(M) + X(Math.min(2 * M, 50000))) / 2, ly + 50, {align: 'center', color: th.time, size: 13, alpha: al});
          K.text(`x = 10${['⁻²', '⁻³', '⁻⁵', '⁻⁹', '⁻¹⁷'][i]}`, (X(M) + X(Math.min(2 * M, 50000))) / 2, ly + 72, {align: 'center', color: th.proved, size: 13, alpha: al});
        });
        K.text('chain length n (log scale): the windows [Mⱼ, 2Mⱼ] cover every n ≥ 1920', lx, ly - 76, {color: th.mute, size: 14, alpha: B});
        K.text('in window j:  thermal purity ≥ 1 − 4xⱼ > ½  ⇒  unique ground state', W / 2, 370, {align: 'center', size: 17, alpha: B * seg(t, 15.5, 16.5)});
        K.tex('\\Delta_n>\\frac{\\log 10}{392}\\approx0.0059', W / 2, 450, {align: 'center', size: 30, color: th.proved, alpha: seg(t, 18.5, 19.5)});
      }
    }
  });

  for (const id of ['ocyl', 'ovec', 'over', 'ofold', 'omoves', 'ostart', 'ogap']) slow(OPEN, id, 1.3);
  root.STBV.OPEN = OPEN;
})(typeof window !== 'undefined' ? window : globalThis);
