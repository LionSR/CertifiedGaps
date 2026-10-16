// Film player shared by film.html and film_open.html. The page sets window.FILM = {list, prefix,
// title, sub} before loading this script; it exposes window.draw(t), window.TOTAL, and the promise
// window.READY that render.mjs waits on before drawing frames.
(function () {
  'use strict';
  const {list: LIST, prefix: PREFIX, title: TITLE, sub: SUB} = window.FILM;
  const {seg, beatAt, texImage, kit} = STBV;
  const S = 2, cv = document.getElementById('c'), ctx = cv.getContext('2d');
  const SERIF = '"EB Garamond", Palatino, serif';
  const TH = {bg: '#fffff8', surface: '#fffff8', ink: '#111111', mute: '#6b6b66', line: '#d9d6c6',
    space: '#a2322e', time: '#34628f', proved: '#8f6412', ok: '#2f6b3d', fail: '#a2322e',
    spaceF: 'rgba(162,50,46,.08)', timeF: 'rgba(52,98,143,.08)', provedF: 'rgba(143,100,18,.10)',
    torusFill: '#efede2'};

  // Every frame must be drawable synchronously, so all TeX is typeset up front: dry-run each
  // chapter to collect its requests, then load the images.
  const TEXC = new Map();
  const texImg = (src, color, px) => TEXC.get(src + '|' + color + '|' + px) || null;
  async function preload() {
    const want = new Map();
    const probe = (src, color, px) => { want.set(src + '|' + color + '|' + px, [src, color, px]); return null; };
    for (const ch of LIST) for (let t = 0; t <= ch.dur; t += .5) { ctx.save(); ch.draw(kit(ctx, TH, probe), t); ctx.restore(); }
    await MathJax.startup.promise;
    await Promise.all([...want].map(async ([key, args]) => TEXC.set(key, await texImage(...args))));
  }

  const INTRO = 5;
  const starts = []; let acc = INTRO; for (const ch of LIST) { starts.push(acc); acc += ch.dur; }
  const TOTAL = acc;

  function caption(text, alpha) {
    ctx.save(); ctx.globalAlpha = alpha; ctx.font = `italic 20px ${SERIF}`; ctx.textAlign = 'center';
    const lines = []; let cur = '';
    for (const w of text.split(' ')) { const tryL = cur ? cur + ' ' + w : w; if (ctx.measureText(tryL).width > 820) { lines.push(cur); cur = w; } else cur = tryL; }
    lines.push(cur);
    const y0 = 487 + (52 - (26 * lines.length + 14)) / 2;
    ctx.fillStyle = TH.ink; lines.forEach((l, i) => ctx.fillText(l, 480, y0 + 26 + i * 26)); ctx.restore();
  }

  function draw(t) {
    ctx.setTransform(S, 0, 0, S, 0, 0); ctx.fillStyle = TH.bg; ctx.fillRect(0, 0, 960, 540);
    if (t < INTRO) {
      const a = Math.min(seg(t, 0, .8), 1 - seg(t, INTRO - .5, INTRO));
      ctx.save(); ctx.globalAlpha = a; ctx.fillStyle = TH.ink; ctx.textAlign = 'center';
      ctx.font = `400 54px ${SERIF}`; ctx.fillText(TITLE, 480, 240);
      ctx.font = `italic 24px ${SERIF}`; ctx.fillText(SUB, 480, 285); ctx.restore();
    } else {
      let i = starts.length - 1; while (i > 0 && t < starts[i]) i--;
      const ch = LIST[i], lt = t - starts[i], fade = Math.min(seg(lt, 0, .4), 1 - seg(lt, ch.dur - .4, ch.dur));
      ctx.save(); ctx.globalAlpha = fade; ctx.translate(72, 18); ctx.scale(.85, .85); ch.draw(kit(ctx, TH, texImg), Math.min(lt, ch.dur)); ctx.restore();
      ctx.save(); ctx.globalAlpha = fade * (1 - seg(lt, 3.5, 4.5) * .55); ctx.fillStyle = TH.mute; ctx.font = `italic 16px ${SERIF}`;
      ctx.fillText(`${PREFIX}${i + 1}  ${ch.title}`, 24, 30); ctx.restore();
      const bi = beatAt(ch, lt), b0 = ch.beats[bi][0], nb = ch.beats[bi + 1] ? ch.beats[bi + 1][0] : ch.dur;
      caption(ch.beats[bi][1], fade * Math.min(seg(lt, b0, b0 + .3), 1 - seg(lt, nb - .25, nb)));
    }
    ctx.fillStyle = TH.ink; ctx.fillRect(0, 539, 960 * t / TOTAL, 1);
  }

  window.TOTAL = TOTAL; window.draw = draw;
  window.READY = (async () => { if (document.fonts) await document.fonts.ready; await preload(); draw(0); })();
})();
