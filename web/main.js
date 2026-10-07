import { loadWasm, WasmCavity, JsCavity } from "./solver.js";

const $ = (id) => document.getElementById(id);
const L = 2; // cavity side length, as in the notebook
const RHO = 1;

// ---- colour maps ---------------------------------------------------------

// Polynomial fit to matplotlib's viridis (the notebook's colormap).
function viridis(t) {
  const c = [
    [0.2777273272234177, 0.005407344544966578, 0.3340998053353061],
    [0.1050930431085774, 1.404613529898575, 1.384590162594685],
    [-0.3308618287255563, 0.214847559468213, 0.09509516302823659],
    [-4.634230498983486, -5.799100973351585, -19.33244095627987],
    [6.228269936347081, 14.17993336680509, 56.69055260068105],
    [4.776384997670288, -13.74514537774601, -65.35303263337234],
    [-5.435455855934631, 4.645852612178535, 26.3124352495832],
  ];
  return [0, 1, 2].map((ch) => {
    let v = 0;
    for (let k = 6; k >= 0; k--) v = v * t + c[k][ch];
    return Math.round(255 * Math.min(1, Math.max(0, v)));
  });
}

// Diverging blue / near-white / red for signed vorticity.
function diverging(t) {
  const stops = [[33, 102, 172], [146, 197, 222], [240, 240, 236], [244, 165, 130], [178, 24, 43]];
  const x = t * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  return stops[i].map((a, ch) => Math.round(a + (stops[i + 1][ch] - a) * f));
}

const lut = (fn) => {
  const out = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) out.set(fn(i / 255), i * 3);
  return out;
};
const MAPS = { viridis: lut(viridis), diverging: lut(diverging) };

// ---- state ---------------------------------------------------------------

const state = {
  wasm: null,
  solver: null,
  engine: "wasm",
  n: 41,
  running: !matchMedia("(prefers-reduced-motion: reduce)").matches,
  time: 0,
  steps: 0,
  msPerStep: 0,
  stepsPerSec: 0,
  view: "p",
};

function params() {
  const s = Number($("nu").value);
  const nu = 0.1 * Math.pow(20, (s - 0.5) * 2); // 0.005 … 2, notebook value at the centre
  const lid = Number($("lid").value);
  const nit = Number($("nit").value);
  const h = L / (state.n - 1);
  // Notebook Δt, reduced when explicit diffusion or advection would go unstable.
  const dt = Math.min(0.001, (0.2 * h * h) / nu, (0.25 * h) / lid);
  return { dt, rho: RHO, nu, lid, nit };
}

function makeSolver(engine, n) {
  return engine === "wasm" && state.wasm ? new WasmCavity(state.wasm, n, n, L, L) : new JsCavity(n, n, L, L);
}

function rebuild({ keepFlow }) {
  const old = state.solver;
  const next = makeSolver(state.engine, state.n);
  next.setParams(params());
  if (keepFlow && old && old.nx === state.n) {
    const a = old.fields();
    const b = next.fields();
    b.u.set(a.u);
    b.v.set(a.v);
    b.p.set(a.p);
  } else {
    state.time = 0;
    state.steps = 0;
  }
  old?.free();
  state.solver = next;
  state.msPerStep = 0;
  seedTracers();
}

// ---- canvases ------------------------------------------------------------

const field = $("field");
const overlay = $("overlay");
const fctx = field.getContext("2d");
const octx = overlay.getContext("2d");
const img = document.createElement("canvas");
const ictx = img.getContext("2d");
let size = 0;

function resize() {
  const css = field.getBoundingClientRect().width;
  const px = Math.max(1, Math.round(css * Math.min(devicePixelRatio || 1, 2)));
  if (px === size) return;
  size = field.width = field.height = overlay.width = overlay.height = px;
  draw();
}

function vorticity(u, v, n) {
  const h = L / (n - 1);
  const w = new Float64Array(n * n);
  for (let j = 1; j < n - 1; j++) {
    for (let i = 1; i < n - 1; i++) {
      const k = j * n + i;
      w[k] = (v[k + 1] - v[k - 1]) / (2 * h) - (u[k + n] - u[k - n]) / (2 * h);
    }
  }
  return w;
}

// Range from the 1st to 99th percentile, so the singular lid corners don't wash out the map.
function robustRange(a) {
  const stride = Math.max(1, Math.floor(a.length / 4000));
  const s = [];
  for (let i = 0; i < a.length; i += stride) s.push(a[i]);
  s.sort((x, y) => x - y);
  return [s[Math.floor(s.length * 0.01)], s[Math.ceil(s.length * 0.99) - 1]];
}

function scalarField() {
  const { u, v, p } = state.solver.fields();
  const n = state.n;
  if (state.view === "speed") {
    const s = new Float64Array(n * n);
    for (let k = 0; k < s.length; k++) s[k] = Math.hypot(u[k], v[k]);
    return { data: s, lo: 0, hi: Math.max(1e-9, Number($("lid").value)), map: "viridis", name: "Speed |u|" };
  }
  if (state.view === "vort") {
    const w = vorticity(u, v, n);
    const [a, b] = robustRange(w);
    const m = Math.max(Math.abs(a), Math.abs(b), 1e-9);
    return { data: w, lo: -m, hi: m, map: "diverging", name: "Vorticity ω = ∂v/∂x − ∂u/∂y" };
  }
  let [lo, hi] = robustRange(p);
  if (hi - lo < 1e-12) hi = lo + 1e-12;
  return { data: p, lo, hi, map: "viridis", name: "Pressure p" };
}

const fmt = (x) => (Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 1 ? x.toFixed(2) : x.toFixed(3));

function draw() {
  if (!state.solver || !size) return;
  const n = state.n;
  const { data, lo, hi, map, name } = scalarField();
  const cmap = MAPS[map];

  if (img.width !== n) img.width = img.height = n;
  const im = ictx.createImageData(n, n);
  const scale = 255 / (hi - lo);
  for (let j = 0; j < n; j++) {
    const row = (n - 1 - j) * n; // y up: lid row at the top of the image
    for (let i = 0; i < n; i++) {
      const c = Math.max(0, Math.min(255, Math.round((data[j * n + i] - lo) * scale))) * 3;
      const o = (row + i) * 4;
      im.data[o] = cmap[c];
      im.data[o + 1] = cmap[c + 1];
      im.data[o + 2] = cmap[c + 2];
      im.data[o + 3] = 255;
    }
  }
  ictx.putImageData(im, 0, 0);
  fctx.imageSmoothingEnabled = true;
  fctx.imageSmoothingQuality = "high";
  // Grid nodes sit on the walls, so map node centres to the canvas edges.
  const cell = size / (n - 1);
  fctx.drawImage(img, -cell / 2, -cell / 2, size + cell, size + cell);

  if ($("arrows").checked) drawArrows();
  drawLegend(lo, hi, map, name);
}

function drawArrows() {
  const { u, v } = state.solver.fields();
  const n = state.n;
  const every = Math.max(1, Math.round((n - 1) / 20));
  const U = Math.max(1e-9, Number($("lid").value));
  const gap = (size / (n - 1)) * every;
  fctx.strokeStyle = "rgba(255,255,255,0.85)";
  fctx.fillStyle = "rgba(255,255,255,0.85)";
  fctx.lineWidth = Math.max(1, size / 500);
  for (let j = every; j < n - 1; j += every) {
    for (let i = every; i < n - 1; i += every) {
      const k = j * n + i;
      const x = (i / (n - 1)) * size;
      const y = size - (j / (n - 1)) * size;
      const len = (Math.hypot(u[k], v[k]) / U) * gap * 0.95;
      if (len < 1) continue;
      const ang = Math.atan2(-v[k], u[k]);
      fctx.save();
      fctx.translate(x, y);
      fctx.rotate(ang);
      fctx.beginPath();
      fctx.moveTo(-len / 2, 0);
      fctx.lineTo(len / 2, 0);
      fctx.stroke();
      const hd = Math.min(len * 0.4, gap * 0.25);
      fctx.beginPath();
      fctx.moveTo(len / 2, 0);
      fctx.lineTo(len / 2 - hd, -hd * 0.5);
      fctx.lineTo(len / 2 - hd, hd * 0.5);
      fctx.fill();
      fctx.restore();
    }
  }
}

const cbar = $("colorbar").getContext("2d");
function drawLegend(lo, hi, map, name) {
  const im = cbar.createImageData(256, 1);
  const cmap = MAPS[map];
  for (let i = 0; i < 256; i++) {
    im.data.set([cmap[i * 3], cmap[i * 3 + 1], cmap[i * 3 + 2], 255], i * 4);
  }
  cbar.putImageData(im, 0, 0);
  $("legendMin").textContent = fmt(lo);
  $("legendMax").textContent = fmt(hi);
  $("legendName").textContent = name;
}

// ---- tracer particles ----------------------------------------------------

const TRACERS = 1400;
const tracers = new Float64Array(TRACERS * 3); // x, y, age

function spawn(i) {
  tracers[i * 3] = 0.02 + Math.random() * (L - 0.04);
  tracers[i * 3 + 1] = 0.02 + Math.random() * (L - 0.04);
  tracers[i * 3 + 2] = Math.floor(Math.random() * 240);
}
function seedTracers() {
  for (let i = 0; i < TRACERS; i++) spawn(i);
  octx.clearRect(0, 0, overlay.width, overlay.height);
}

function sample(a, x, y) {
  const n = state.n;
  const h = L / (n - 1);
  const fx = Math.min(n - 1.0001, Math.max(0, x / h));
  const fy = Math.min(n - 1.0001, Math.max(0, y / h));
  const i = Math.floor(fx), j = Math.floor(fy);
  const tx = fx - i, ty = fy - j;
  const k = j * n + i;
  return (a[k] * (1 - tx) + a[k + 1] * tx) * (1 - ty) + (a[k + n] * (1 - tx) + a[k + n + 1] * tx) * ty;
}

function moveTracers(dtSim) {
  const { u, v } = state.solver.fields();
  const px = size / L;
  octx.globalCompositeOperation = "destination-out";
  octx.fillStyle = "rgba(0,0,0,0.12)";
  octx.fillRect(0, 0, size, size);
  octx.globalCompositeOperation = "source-over";
  if (!$("tracers").checked) return;
  octx.strokeStyle = "rgba(255,255,255,0.8)";
  octx.lineWidth = Math.max(1, size / 420);
  octx.beginPath();
  for (let i = 0; i < TRACERS; i++) {
    const o = i * 3;
    const x = tracers[o], y = tracers[o + 1];
    // Midpoint (RK2) step through the bilinear velocity field.
    const xm = x + 0.5 * dtSim * sample(u, x, y);
    const ym = y + 0.5 * dtSim * sample(v, x, y);
    const nx = x + dtSim * sample(u, xm, ym);
    const ny = y + dtSim * sample(v, xm, ym);
    tracers[o + 2] -= 1;
    if (tracers[o + 2] < 0 || nx <= 0 || nx >= L || ny <= 0 || ny >= L) {
      spawn(i);
      continue;
    }
    octx.moveTo(x * px, size - y * px);
    octx.lineTo(nx * px, size - ny * px);
    tracers[o] = nx;
    tracers[o + 1] = ny;
  }
  octx.stroke();
}

// ---- main loop -----------------------------------------------------------

const FRAME_BUDGET_MS = 14;

function advance() {
  const s = state.solver;
  const target = Number($("spf").value);
  const t0 = performance.now();
  let done = 0;
  // Run in small chunks so a slow engine or big grid can't freeze the page.
  while (done < target) {
    const chunk = Math.min(target - done, Math.max(1, Math.ceil(target / 6)));
    s.step(chunk);
    done += chunk;
    if (performance.now() - t0 > FRAME_BUDGET_MS) break;
  }
  const ms = performance.now() - t0;
  const per = ms / done;
  state.msPerStep = state.msPerStep ? state.msPerStep * 0.9 + per * 0.1 : per;
  state.steps += done;
  const dt = params().dt;
  state.time += done * dt;
  return done * dt;
}

function updateReadout() {
  const p = params();
  $("re").textContent = ((p.lid * L) / p.nu).toFixed(p.lid * L / p.nu < 10 ? 1 : 0);
  $("time").textContent = state.time.toFixed(3);
  $("dt").textContent = p.dt < 0.001 ? p.dt.toExponential(2) : p.dt.toFixed(3);
  $("steps").textContent = state.steps.toLocaleString("en");
  $("ms").textContent = state.msPerStep ? state.msPerStep.toFixed(3) : "–";
  $("sps").textContent = state.msPerStep ? Math.round(1000 / state.msPerStep).toLocaleString("en") : "–";
}

function frame() {
  if (state.running) {
    const dtSim = advance();
    moveTracers(dtSim);
  }
  draw();
  updateReadout();
  requestAnimationFrame(frame);
}

// ---- controls ------------------------------------------------------------

function syncLabels() {
  const p = params();
  $("nuOut").textContent = p.nu.toPrecision(2);
  $("lidOut").textContent = p.lid.toFixed(1);
  $("lidLabel").textContent = p.lid.toFixed(1);
  $("nitOut").textContent = p.nit;
  $("spfOut").textContent = $("spf").value;
}

function onParams() {
  syncLabels();
  state.solver.setParams(params());
}
for (const id of ["nu", "lid", "nit"]) $(id).addEventListener("input", onParams);
$("spf").addEventListener("input", syncLabels);

$("grid").addEventListener("change", () => {
  state.n = Number($("grid").value);
  rebuild({ keepFlow: false });
  syncLabels();
});

for (const r of document.querySelectorAll('input[name="view"]')) {
  r.addEventListener("change", () => {
    state.view = r.value;
  });
}
for (const r of document.querySelectorAll('input[name="engine"]')) {
  r.addEventListener("change", () => {
    state.engine = r.value;
    rebuild({ keepFlow: true });
  });
}
$("tracers").addEventListener("change", () => octx.clearRect(0, 0, size, size));

function setRunning(on) {
  state.running = on;
  $("play").textContent = on ? "Pause" : "Run";
  document.documentElement.style.setProperty("--lid-play", on ? "running" : "paused");
}
$("play").addEventListener("click", () => setRunning(!state.running));
$("reset").addEventListener("click", () => {
  state.solver.reset();
  state.time = 0;
  state.steps = 0;
  seedTracers();
});

// ---- start ---------------------------------------------------------------

async function start() {
  try {
    const bytes = await (await fetch(new URL("cavity.wasm", import.meta.url))).arrayBuffer();
    state.wasm = await loadWasm(bytes);
  } catch (err) {
    console.warn("WebAssembly unavailable, using the JavaScript solver", err);
    state.engine = "js";
    $("eJ").checked = true;
    $("eW").disabled = true;
    $("engineNote").textContent = "This browser couldn't load the WebAssembly module, so the JavaScript solver is running.";
  }
  syncLabels();
  rebuild({ keepFlow: false });
  // Open on a developed flow rather than fluid at rest.
  state.solver.step(400);
  state.steps = 400;
  state.time = 400 * params().dt;
  setRunning(state.running);
  new ResizeObserver(resize).observe(field);
  resize();
  requestAnimationFrame(frame);
}
start();
