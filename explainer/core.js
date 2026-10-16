// Shared core for the spacetime-bootstrap page and film.
// (1) exact two-purity recurrence in log space (mirrors julia/src/recurrence/map.jl),
// (2) cancellation-free TFIM purity defects from free fermions,
// (3) a small painter's-algorithm renderer for the spacetime sheet -> torus morph.
(function (root) {
  'use strict';
  const LN2 = Math.log(2), ZSTAR = 2 - Math.SQRT2;
  const lae = (a, b) => { const m = Math.max(a, b); return m === -Infinity ? m : m + Math.log(Math.exp(a - m) + Math.exp(b - m)); };

  // ---------- recurrence ----------
  // z(a,b) = 1 - (1-b)(1-a)^2 from log-defects, as [z, log z]
  function zlog(la, lb) {
    if (la < -30 && lb < -30) { const lz = lae(lb, LN2 + la); return [Math.exp(lz), lz]; }
    const z = -Math.expm1(Math.log1p(-Math.exp(lb)) + 2 * Math.log1p(-Math.exp(la)));
    return [z, Math.log(z)];
  }
  const fLog = (z, lz) => 2 * lz - LN2 - 2 * Math.log1p(-z);
  // Orbit of the exact map: rows {j, ls, lt, z, lz, rate}; rate is beta0 * b_j.
  function orbit(es, et, steps = 60) {
    const rows = [];
    if (!(es > 0 && et > 0)) return {rows, ok: false};
    let ls = Math.log(es), lt = Math.log(et);
    for (let j = 0; j < steps; j++) {
      const [z, lz] = zlog(ls, lt);
      rows.push({j, ls, lt, z, lz, rate: z < .5 ? (Math.log1p(-z) - lz) / Math.pow(2, j) : -Infinity});
      if (z >= ZSTAR) return {rows, ok: false};
      const ls2 = fLog(z, lz);
      const [z2, lz2] = zlog(lt, ls2);
      if (z2 >= ZSTAR) return {rows, ok: false};
      lt = fLog(z2, lz2); ls = ls2;
    }
    return {rows, ok: true};
  }
  // beta0 * (uniform gap bound), or null if the seed does not close the recurrence
  function scaledBound(es, et) {
    const {rows, ok} = orbit(es, et);
    if (!ok) return null;
    let m = Infinity; for (const r of rows) m = Math.min(m, r.rate);
    return m > 0 && isFinite(m) ? m : null;
  }

  // ---------- TFIM H = -sum Z_j Z_{j+1} - g sum X_j on a ring ----------
  const disp = (k, g) => 2 * Math.sqrt(1 + g * g - 2 * g * Math.cos(k));
  // Fourier coefficient (1/2pi) int disp(k) cos(qk) dk, from disp = 2 max(1,g) |1 - c e^{ik}|
  function dispHat(q, g) {
    const c = Math.min(g, 1 / g), pref = 2 * Math.max(1, g), c2 = c * c;
    let aq = 1; for (let m = 0; m < q; m++) aq *= (m - 0.5) / (m + 1);
    let am = 1, amq = aq, s = 0, w = 1;
    for (let m = 0; m < 200000; m++) {
      const term = am * amq * w; s += term;
      if (m > 10 && Math.abs(term) < 1e-18 * Math.abs(s)) break;
      am *= (m - 0.5) / (m + 1); amq *= (m + q - 0.5) / (m + q + 1); w *= c2;
    }
    return pref * s * Math.pow(c, q);
  }
  // E_NS(2n) - 2 E_NS(n), exponentially small: Poisson summation over the NS momenta
  function mismatch(n, g) {
    let d = 0;
    for (let r = 1; r < 4000; r++) {
      const t = dispHat(r * n, g);
      if (t === 0 || (r > 2 && Math.abs(t) < 1e-19 * Math.abs(d))) break;
      d += 4 * (r % 2 ? -1 : 1) * (dispHat(2 * r * n, g) - t);
    }
    return -n * d / 2;
  }
  // u = Z_n(beta) e^{beta E_NS(n)} - 1, formed without cancellation
  function uExc(n, beta, g) {
    let a = 0, b = 0, ER = 0, ENS = 0, cR = 0, dR = 0, sgn = 1;
    for (let m = 0; m < n; m++) {
      const eN = disp(2 * Math.PI * (m + .5) / n, g), qN = Math.exp(-beta * eN);
      ENS -= eN / 2; a += Math.log1p(qN); b += Math.log1p(-qN);
      const eR = m === 0 ? 2 * (g - 1) : disp(2 * Math.PI * m / n, g);
      if (eR < 0) sgn = -sgn;
      const aR = Math.abs(eR), qR = Math.exp(-beta * aR);
      ER -= aR / 2; cR += Math.log1p(qR); dR += aR === 0 ? -Infinity : Math.log1p(-qR);
    }
    const ns = 0.5 * (Math.expm1(a) + Math.expm1(b));
    const w = Math.exp(-beta * (ER - ENS));
    const rr = sgn > 0
      ? 0.5 * w * (dR === -Infinity ? Math.exp(cR) : Math.exp(dR) * Math.expm1(cR - dR))
      : 0.5 * w * (Math.exp(cR) + Math.exp(dR));
    return ns + rr;
  }
  // [1 - P_s(n, beta), 1 - P_tau(n, beta)]
  // Each defect is floored at a generous estimate of its rounding error, so a float-noise
  // value never pretends to be a tiny defect.
  function tfimDefects(n, beta, g) {
    const u1 = Math.log1p(uExc(n, beta, g)), u2 = Math.log1p(uExc(2 * n, beta, g)), u3 = Math.log1p(uExc(n, 2 * beta, g));
    const ls = -beta * mismatch(n, g) + u2 - 2 * u1, lt = u3 - 2 * u1;
    const tol = x => 64 * Number.EPSILON * x + 1e-300;
    return [Math.max(-Math.expm1(ls), tol(u2 + 2 * u1)), Math.max(-Math.expm1(lt), tol(u3 + 2 * u1))];
  }
  // Seed defects of the criterion: eps_s = 1 - P_s(n0, b0), eps_t = 1 - P_tau(2 n0, b0)
  function tfimSeed(n0, b0, g) {
    return [tfimDefects(n0, b0, g)[0], tfimDefects(2 * n0, b0, g)[1]];
  }
  const tfimGap = g => 2 * Math.abs(1 - g);           // infinite-chain gap
  const tfimXi = g => 1 / Math.abs(Math.log(g));      // correlation length

  // ---------- spacetime sheet -> torus ----------
  // m in [0,2]: 0 flat sheet, 1 cylinder (time rolled), 2 torus (space rolled too).
  function sheetPoint(u, v, m, Wd, Hd) {
    const m1 = Math.min(1, Math.max(0, m)), m2 = Math.min(1, Math.max(0, m - 1));
    const k1 = Math.max(1e-6, m1) * 2 * Math.PI / Hd, k2 = Math.max(1e-6, m2) * 2 * Math.PI / Wd;
    const y = Math.sin(k1 * v) / k1, z = (1 - Math.cos(k1 * v)) / k1, zc = Hd / (2 * Math.PI);
    const R = 1 / k2 + zc - z;
    // recentre: flat sheet sits at z = -zc; a partial bend of angle th moves the arc centroid
    const th = k2 * Wd / 2, arc = (1 - (th < 1e-9 ? 1 : Math.sin(th) / th)) / k2;
    return [Math.sin(k2 * u) * R, y, -(Math.cos(k2 * u) * R - 1 / k2) + zc * (1 - m1) - arc];
  }
  // Draw the sheet. opt: {cx, cy, scale, m, yaw, pitch, cols, rows, colHi (column index or -1),
  //   colors: {space, time, fill, hi, edge}, lineW}
  function drawTorus(ctx, opt) {
    const Wd = 2, Hd = 1, NU = opt.nu || 64, NV = opt.nv || 32;
    const cy = Math.cos(opt.yaw), sy = Math.sin(opt.yaw), cp = Math.cos(opt.pitch), sp = Math.sin(opt.pitch);
    const P = (u, v) => {
      let [x, y, z] = sheetPoint(u, v, opt.m, Wd, Hd);
      const x1 = cy * x + sy * z, z1 = -sy * x + cy * z;
      const y2 = cp * y - sp * z1, z2 = sp * y + cp * z1;
      return [x1, y2, z2];
    };
    const pts = [];
    for (let i = 0; i <= NU; i++) { pts.push([]); for (let j = 0; j <= NV; j++) pts[i].push(P(-Wd / 2 + Wd * i / NU, -Hd / 2 + Hd * j / NV)); }
    const quads = [];
    for (let i = 0; i < NU; i++) for (let j = 0; j < NV; j++) {
      const a = pts[i][j], b = pts[i + 1][j], c = pts[i + 1][j + 1], d = pts[i][j + 1];
      const ux = c[0] - a[0], uy = c[1] - a[1], uz = c[2] - a[2], vx = d[0] - b[0], vy = d[1] - b[1], vz = d[2] - b[2];
      let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      const nl = Math.hypot(nx, ny, nz) || 1; nx /= nl; ny /= nl; nz /= nl;
      quads.push({i, j, q: [a, b, c, d], depth: (a[2] + b[2] + c[2] + d[2]) / 4, shade: Math.abs(-0.35 * nx - 0.55 * ny + 0.76 * nz)});
    }
    quads.sort((p, q) => p.depth - q.depth);
    const X = p => opt.cx + p[0] * opt.scale, Y = p => opt.cy - p[1] * opt.scale;
    const cols = opt.cols || 16, rows = opt.rows || 8, col = opt.colors;
    const hiCol = opt.colHi ?? -1;
    ctx.save(); ctx.lineJoin = 'round';
    for (const Q of quads) {
      const inHi = hiCol >= 0 && Math.floor(Q.i / (NU / cols)) === hiCol;
      ctx.beginPath(); ctx.moveTo(X(Q.q[0]), Y(Q.q[0]));
      for (let k = 1; k < 4; k++) ctx.lineTo(X(Q.q[k]), Y(Q.q[k]));
      ctx.closePath();
      ctx.globalAlpha = 1;
      ctx.fillStyle = inHi ? col.hi : col.fill;
      ctx.fill();
      ctx.fillStyle = `rgba(0,0,0,${(1 - Q.shade) * (opt.shadeAmt ?? .45)})`; ctx.fill();
      // grid lines on this quad's edges
      ctx.lineWidth = opt.lineW || 1.2;
      if (Q.j % (NV / rows) === 0) { ctx.strokeStyle = col.space; ctx.beginPath(); ctx.moveTo(X(Q.q[0]), Y(Q.q[0])); ctx.lineTo(X(Q.q[1]), Y(Q.q[1])); ctx.stroke(); }
      if (Q.i % (NU / cols) === 0) { ctx.strokeStyle = col.time; ctx.beginPath(); ctx.moveTo(X(Q.q[0]), Y(Q.q[0])); ctx.lineTo(X(Q.q[3]), Y(Q.q[3])); ctx.stroke(); }
      if (Q.j === NV - 1) { ctx.strokeStyle = col.space; ctx.beginPath(); ctx.moveTo(X(Q.q[3]), Y(Q.q[3])); ctx.lineTo(X(Q.q[2]), Y(Q.q[2])); ctx.stroke(); }
      if (Q.i === NU - 1) { ctx.strokeStyle = col.time; ctx.beginPath(); ctx.moveTo(X(Q.q[1]), Y(Q.q[1])); ctx.lineTo(X(Q.q[2]), Y(Q.q[2])); ctx.stroke(); }
    }
    ctx.restore();
  }

  root.STB = {orbit, scaledBound, tfimDefects, tfimSeed, tfimGap, tfimXi, mismatch, drawTorus, sheetPoint, ZSTAR};
})(typeof window !== 'undefined' ? window : globalThis);
