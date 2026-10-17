// Render one of the film pages to an MP4 with a spoken narration track:
//   node render.mjs film.html film.mp4 [--mute]
// Needs Playwright (with Chromium), kokoro-js (unless --mute), and ffmpeg on the PATH. Frames are
// drawn deterministically by window.draw(t), so the output does not depend on machine speed.
//
// Each caption beat is read aloud from narration.js by the Kokoro-82M speech model, run locally.
// When a line outlasts its beat, the beat plays up to 25% slower and then holds its last full
// frame until the line ends, so the picture never runs ahead of the voice. Clips are cached in
// .voice/ by a hash of their text, so editing one line re-synthesizes only that line.
import { chromium } from 'playwright';
import { spawn } from 'child_process';
import { createHash } from 'crypto';
import fs from 'fs';
import path from 'path';

const args = process.argv.slice(2), MUTE = args.includes('--mute');
const [page = 'film.html', out = 'film.mp4'] = args.filter(a => !a.startsWith('--'));
const FPS = 30, SR = 24000, MODEL = 'onnx-community/Kokoro-82M-v1.0-ONNX', VOICE = 'af_heart', SPEED = 1;
const LEAD = .35, TRAIL = .55, TAIL = .5, MAXSLOW = 1.25;  // seconds; TAIL keeps fades out of a hold

const browser = await chromium.launch();
const tab = await browser.newPage({ viewport: { width: 1920, height: 1080 } });
tab.on('pageerror', e => console.error('page error:', e.message));
await tab.goto('file://' + path.resolve(page));
await tab.evaluate(() => window.READY);  // set synchronously by film.js; resolves once TeX is loaded
const total = await tab.evaluate(() => window.TOTAL);
const beats = await tab.evaluate(() => window.BEATS);

// ---- narration clips (Float32 mono at SR), cached on disk ----
const CACHE = path.resolve(path.dirname(page), '.voice');
let tts = null;
// Kokoro silently truncates long inputs, so a long line is voiced in sentence-sized chunks of at
// most CHUNK characters, joined by short pauses. Short lines keep their original cache keys.
const CHUNK = 300, GAP = .22;
function chunks(text) {
  if (text.length <= CHUNK) return [text];
  const out = []; let cur = '';
  for (const s of text.split(/(?<=[.?!:;])\s+/)) { if (cur && cur.length + s.length + 1 > CHUNK) { out.push(cur); cur = s; } else cur = cur ? cur + ' ' + s : s; }
  return [...out, cur];
}
async function clip(text) {
  const parts = chunks(text), key = createHash('sha256').update(JSON.stringify(parts.length > 1 ? [MODEL, VOICE, SPEED, parts] : [MODEL, VOICE, SPEED, text])).digest('hex').slice(0, 24);
  const file = path.join(CACHE, key + '.f32');
  if (!fs.existsSync(file)) {
    if (!tts) tts = await (await import('kokoro-js')).KokoroTTS.from_pretrained(MODEL, { dtype: 'fp32', device: 'cpu' });
    const pieces = [];
    for (const p of parts) {
      const a = await tts.generate(p, { voice: VOICE, speed: SPEED });
      if (a.sampling_rate !== SR) throw new Error(`unexpected sample rate ${a.sampling_rate}`);
      if (pieces.length) pieces.push(new Float32Array(Math.round(GAP * SR)));
      pieces.push(a.audio);
    }
    const all = new Float32Array(pieces.reduce((n, p) => n + p.length, 0)); let o = 0;
    for (const p of pieces) { all.set(p, o); o += p.length; }
    fs.mkdirSync(CACHE, { recursive: true });
    fs.writeFileSync(file, Buffer.from(all.buffer));
  }
  const b = fs.readFileSync(file);
  return new Float32Array(b.buffer, b.byteOffset, b.byteLength / 4);
}
for (const b of beats) {
  if (b.mismatch) throw new Error(`narration.js: ${b.where.split('[')[0]} has the wrong number of lines`);
  if (!MUTE && typeof b.say !== 'string') throw new Error(`narration.js: no line for ${b.where || 'the title card'}`);
}
const clips = [];
for (const [i, b] of beats.entries()) {
  clips.push(MUTE ? new Float32Array(0) : await clip(b.say));
  process.stdout.write(`\rvoice ${i + 1}/${beats.length}`);
}
console.log();

// ---- schedule: output duration D of each beat and the output-to-film time map ----
// A beat of film length L is shown for D >= L seconds: its first L - h seconds play slowed by s,
// then the frame at t1 - h holds for F seconds, then the last h seconds play (fades included).
let o = 0;
const sched = beats.map((b, i) => {
  const L = b.t1 - b.t0, h = Math.min(TAIL, L / 3), need = clips[i].length / SR;
  const D = need ? Math.max(L, LEAD + need + TRAIL + h) : L, s = Math.min(D / L, MAXSLOW);
  const r = { o0: o, D, f0: b.t0, f1: b.t1, L, h, s, F: D - L * s };
  o += D; return r;
});
const OUT = o;
function filmTime(u) {
  let k = 0; while (k < sched.length - 1 && u >= sched[k].o0 + sched[k].D) k++;
  const { o0, f0, f1, L, h, s, F } = sched[k], x = u - o0, A = (L - h) * s;
  return x < A ? f0 + x / s : x < A + F ? f1 - h : Math.min(f1, f1 - h + (x - A - F) / s);
}

// ---- audio: place each clip LEAD seconds into its beat, then write a 16-bit WAV ----
const track = new Float32Array(Math.ceil(OUT * SR));
sched.forEach((r, i) => track.set(clips[i].subarray(0, track.length - Math.round((r.o0 + LEAD) * SR)), Math.round((r.o0 + LEAD) * SR)));
const wav = Buffer.alloc(44 + track.length * 2);
wav.write('RIFF', 0); wav.writeUInt32LE(36 + track.length * 2, 4); wav.write('WAVEfmt ', 8);
wav.writeUInt32LE(16, 16); wav.writeUInt16LE(1, 20); wav.writeUInt16LE(1, 22); wav.writeUInt32LE(SR, 24);
wav.writeUInt32LE(SR * 2, 28); wav.writeUInt16LE(2, 32); wav.writeUInt16LE(16, 34); wav.write('data', 36);
wav.writeUInt32LE(track.length * 2, 40);
track.forEach((v, j) => wav.writeInt16LE(Math.round(Math.max(-1, Math.min(1, v)) * 32767), 44 + 2 * j));
const wavFile = out + '.narration.wav';
fs.writeFileSync(wavFile, wav);

// ---- video: frames piped to ffmpeg, muxed with the loudness-normalized narration ----
const audioArgs = MUTE ? ['-an'] : ['-map', '0:v', '-map', '1:a', '-af', 'loudnorm=I=-16:TP=-1.5:LRA=11',
  '-c:a', 'aac', '-b:a', '160k', '-ar', '48000'];
const ff = spawn('ffmpeg', ['-y', '-f', 'image2pipe', '-framerate', String(FPS), '-i', '-', ...(MUTE ? [] : ['-i', wavFile]),
  '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '22', '-preset', 'slow', ...audioArgs, '-movflags', '+faststart', out],
  { stdio: ['pipe', 'ignore', 'inherit'] });
const n = Math.round(OUT * FPS);
let lastT = NaN, png = null;
for (let i = 0; i < n; i++) {
  const u = i / FPS, t = filmTime(u);
  // Held frames differ only in the progress bar, a one-pixel line; reuse them to save time.
  if (t !== lastT || i % FPS === 0) {
    png = Buffer.from(await tab.evaluate(([t, p]) => { draw(t, p); return document.getElementById('c').toDataURL('image/png').split(',')[1]; }, [t, u / OUT]), 'base64');
    lastT = t;
  }
  if (!ff.stdin.write(png)) await new Promise(r => ff.stdin.once('drain', r));
  if (i % (10 * FPS) === 0) process.stdout.write(`\rframes ${i}/${n}`);
}
ff.stdin.end();
await new Promise(r => ff.on('close', r));
await browser.close();
fs.rmSync(wavFile);
console.log(`\r${out}: ${OUT.toFixed(1)} s (film clock ${total.toFixed(1)} s)`);
