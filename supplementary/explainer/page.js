// Shared by index.html and open.html: drawing helpers, theme colors, TeX on canvas, and the
// figure players. Load it after the article; each page then adds its own widgets and pushes
// their redraw functions onto THEMED, so a theme switch repaints them too.
const NS = 'http://www.w3.org/2000/svg';
const $ = id => document.getElementById(id);
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e); return e;
}
function txt(parent, x, y, s, attrs = {}) {
  const t = el('text', Object.assign({x, y, 'font-family': 'IBM Plex Mono, monospace', 'font-size': 12, fill: 'var(--mute)'}, attrs), parent);
  t.textContent = s; return t;
}
function barsSvg(svg, x, y, w, h, vals, color) {
  const gap = 6, bw = (w - gap*(vals.length - 1)) / vals.length;
  el('line', {x1: x, y1: y+h, x2: x+w, y2: y+h, stroke: 'var(--line)'}, svg);
  vals.forEach((v, i) => {
    el('rect', {x: x + i*(bw+gap), y: y + h - v*h, width: bw, height: Math.max(0, v*h), fill: color, opacity: i ? .5 : 1, rx: 2}, svg);
    if (v > .04) txt(svg, x + i*(bw+gap) + bw/2, y + h - v*h - 5, v.toFixed(2), {'text-anchor': 'middle', fill: 'var(--ink)', 'font-size': 11});
  });
}
const sup = n => String(n).split('').map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(c)] || c).join('');
function sci(x, d = 2) {
  if (!isFinite(x)) return '–'; if (x === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(x))), m = x / Math.pow(10, e);
  return Math.abs(e) < 3 ? x.toPrecision(d + 1) : `${m.toFixed(d)}×10${sup(e)}`;
}
const fmtLog = l => { if (!isFinite(l)) return '0'; const e = l / Math.LN10, ei = Math.floor(e); return e > -300 ? sci(Math.exp(l)) : `${Math.pow(10, e - ei).toFixed(2)}×10${sup(ei)}`; };

const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;

// ---------- theme ----------
let TH = null;
function readTheme() {
  const g = v => css(v);
  TH = {bg: g('--stage'), surface: g('--surface'), ink: g('--ink'), mute: g('--mute'), line: g('--line'), space: g('--space'), time: g('--time'),
    proved: g('--proved'), ok: g('--ok'), fail: g('--fail'), spaceF: g('--space-f'), timeF: g('--time-f'), provedF: g('--proved-f'), torusFill: g('--torus-fill')};
}
readTheme();

// ---------- TeX on canvas ----------
const TEXC = new Map(); let texReady = false;
function texImg(src, color, px) {
  const key = src + '|' + color + '|' + px; const hit = TEXC.get(key);
  if (hit) return hit === 'pending' ? null : hit;
  if (!texReady) return null;
  TEXC.set(key, 'pending');
  STBV.texImage(src, color, px).then(im => { TEXC.set(key, im); redrawAll(); }, () => TEXC.delete(key));
  return null;
}

// ---------- players ----------
const PLAYERS = [];
function makePlayer(host) {
  const ch = [...STBV.CH, ...(STBV.OPEN || [])].find(c => c.id === host.dataset.ch); if (!ch) return;
  host.innerHTML = `<canvas width="960" height="540" aria-label="${ch.title} animation"></canvas>
    <div class="r pmeta"><p class="main pcap" aria-live="polite"></p><div class="side"><p class="pctl"><button type="button" class="pp">play</button> · <button type="button" class="pv">back</button> · <button type="button" class="nx">next</button>
    <span class="bt"></span></p><input type="range" min="0" max="${ch.dur}" step="0.01" aria-label="Scrub"></div></div>`;
  const cv = host.querySelector('canvas'), cap = host.querySelector('.pcap'), pp = host.querySelector('.pp'), sl = host.querySelector('input'), bt = host.querySelector('.bt');
  const P = {ch, t: ch.dur, playing: false, seen: false};
  function size() { const dpr = Math.min(2, window.devicePixelRatio || 1); cv.width = 960 * dpr; cv.height = 540 * dpr; P.dpr = dpr; }
  size();
  const beatIdx = () => STBV.beatAt(ch, P.t);
  P.draw = () => {
    const ctx = cv.getContext('2d'); ctx.setTransform(P.dpr, 0, 0, P.dpr, 0, 0);
    ctx.fillStyle = TH.bg; ctx.fillRect(0, 0, 960, 540);
    ch.draw(STBV.kit(ctx, TH, texImg), P.t);
    const bi = beatIdx();
    if (cap.dataset.bi !== String(bi)) { cap.textContent = ch.beats[bi][1]; cap.dataset.bi = bi; }
    sl.value = P.t; bt.textContent = `step ${bi + 1} of ${ch.beats.length}`;
    pp.textContent = P.playing ? 'pause' : (P.t >= ch.dur - .01 ? 'replay' : 'play'); pp.setAttribute('aria-label', P.playing ? 'Pause' : 'Play');
  };
  pp.addEventListener('click', () => { if (P.playing) P.playing = false; else { if (P.t >= ch.dur - .01) P.t = 0; P.playing = true; } P.draw(); });
  host.querySelector('.nx').addEventListener('click', () => { const bi = beatIdx(); P.playing = false; P.t = bi + 1 < ch.beats.length ? ch.beats[bi + 1][0] + 2.5 : ch.dur; P.t = Math.min(P.t, bi + 2 < ch.beats.length ? ch.beats[bi + 2][0] - .01 : ch.dur); P.draw(); });
  host.querySelector('.pv').addEventListener('click', () => { const bi = beatIdx(); P.playing = false; P.t = bi > 0 ? Math.min(ch.beats[bi - 1][0] + 2.5, ch.beats[bi][0] - .01) : 0; P.draw(); });
  sl.addEventListener('input', () => { P.playing = false; P.t = +sl.value; P.draw(); });
  PLAYERS.push(P);
  new IntersectionObserver(es => es.forEach(en => { if (en.isIntersecting && !P.seen && !reduce) { P.seen = true; P.t = 0; P.playing = true; } }), {threshold: .55}).observe(host);
}
document.querySelectorAll('.player').forEach(makePlayer);
function redrawAll() { PLAYERS.forEach(P => P.draw()); }
let last = performance.now();
function tick(now) {
  const dt = Math.min(.05, (now - last) / 1000); last = now;
  for (const P of PLAYERS) if (P.playing) { P.t = Math.min(P.ch.dur, P.t + dt); if (P.t >= P.ch.dur) P.playing = false; P.draw(); }
  requestAnimationFrame(tick);
}
redrawAll(); requestAnimationFrame(tick);
(function waitMJ() { if (window.MathJax && MathJax.startup && MathJax.startup.promise) MathJax.startup.promise.then(() => { texReady = true; redrawAll(); }); else setTimeout(waitMJ, 100); })();
if (document.fonts && document.fonts.ready) document.fonts.ready.then(redrawAll);


// ---------- nav ----------
const links = [...document.querySelectorAll('#toc a')];
const io = new IntersectionObserver(es => es.forEach(en => { if (en.isIntersecting) links.forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + en.target.id)); }), {rootMargin: '-40% 0px -55% 0px'});
document.querySelectorAll('section.ch').forEach(s => io.observe(s));
const THEMED = [];
function themeRedraw() { readTheme(); TEXC.clear(); redrawAll(); THEMED.forEach(f => f()); }
try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', themeRedraw); } catch (e) {}
new MutationObserver(themeRedraw).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
