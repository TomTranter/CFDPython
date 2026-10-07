// Compare the wasm and JavaScript solvers against the notebook's NumPy code.
// Usage (from web/): node validate/check.mjs
import { execFileSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { loadWasm, WasmCavity, JsCavity } from "../solver.js";

const here = (f) => fileURLToPath(new URL(f, import.meta.url));
const wasm = await loadWasm(readFileSync(here("../cavity.wasm")));
const params = { dt: 0.001, rho: 1, nu: 0.1, lid: 1, nit: 50 };
const TOL = 1e-9;
let failed = false;

for (const nt of [100, 700]) {
  const ref = JSON.parse(execFileSync("python3", [here("reference.py"), String(nt)], { maxBuffer: 1 << 26 }));
  const solvers = { wasm: new WasmCavity(wasm, ref.nx, ref.ny), js: new JsCavity(ref.nx, ref.ny) };
  for (const [name, s] of Object.entries(solvers)) {
    s.setParams(params);
    s.step(nt);
    const f = s.fields();
    const errs = ["u", "v", "p"].map((k) => {
      let m = 0;
      for (let i = 0; i < f[k].length; i++) m = Math.max(m, Math.abs(f[k][i] - ref[k][i]));
      return `${k} ${m.toExponential(2)}`;
    });
    const worst = Math.max(...errs.map((e) => Number(e.split(" ")[1])));
    const ok = worst < TOL;
    failed ||= !ok;
    console.log(`nt=${nt} ${name.padEnd(4)} max |diff| vs NumPy: ${errs.join(", ")}  ${ok ? "OK" : "FAIL"}`);
    s.free();
  }
}

// Rough speed comparison on a larger grid.
for (const n of [41, 129]) {
  for (const [name, s] of [["wasm", new WasmCavity(wasm, n, n)], ["js", new JsCavity(n, n)]]) {
    s.setParams(params);
    s.step(5); // warm up
    const t0 = performance.now();
    s.step(100);
    console.log(`${n}x${n} ${name.padEnd(4)} ${((performance.now() - t0) / 100).toFixed(3)} ms/step`);
    s.free();
  }
}
process.exit(failed ? 1 : 0);
