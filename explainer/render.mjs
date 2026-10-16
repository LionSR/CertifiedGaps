// Render one of the film pages to an MP4: node render.mjs film.html film.mp4
// Needs Playwright (with Chromium) and ffmpeg on the PATH. Frames are drawn deterministically by
// window.draw(t), so the output does not depend on machine speed.
import { chromium } from 'playwright';
import { spawn } from 'child_process';
import path from 'path';

const [page = 'film.html', out = 'film.mp4', fps = '30'] = process.argv.slice(2);
const browser = await chromium.launch();
const tab = await browser.newPage({ viewport: { width: 1920, height: 1080 } });
tab.on('pageerror', e => console.error('page error:', e.message));
await tab.goto('file://' + path.resolve(page));
await tab.evaluate(() => window.READY);  // set synchronously by film.js; resolves once TeX is loaded
const total = await tab.evaluate(() => window.TOTAL);
const ff = spawn('ffmpeg', ['-y', '-f', 'image2pipe', '-framerate', fps, '-i', '-', '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
  '-crf', '22', '-preset', 'slow', '-movflags', '+faststart', out], { stdio: ['pipe', 'ignore', 'inherit'] });
const n = Math.round(total * +fps);
for (let i = 0; i < n; i++) {
  const png = await tab.evaluate(t => { draw(t); return document.getElementById('c').toDataURL('image/png').split(',')[1]; }, i / +fps);
  if (!ff.stdin.write(Buffer.from(png, 'base64'))) await new Promise(r => ff.stdin.once('drain', r));
}
ff.stdin.end();
await new Promise(r => ff.on('close', r));
await browser.close();
console.log(`${out}: ${total.toFixed(1)} s`);
