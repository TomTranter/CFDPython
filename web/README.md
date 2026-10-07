# Lid-driven cavity in the browser

An interactive port of [Step 11](../lessons/14_Step_11.ipynb), the lid-driven cavity flow. The finite-difference scheme is the notebook's own, rewritten in Rust and compiled to WebAssembly. A JavaScript copy of the same scheme runs alongside it for speed comparison, and the page falls back to it if WebAssembly is unavailable.

## Run it

The compiled module `cavity.wasm` is committed, so you only need a static file server:

```sh
cd web
python3 -m http.server 8000
# open http://localhost:8000
```

Opening `index.html` directly from disk won't work, because browsers block `fetch` of local files.

## Rebuild the WebAssembly module

```sh
rustup target add wasm32-unknown-unknown   # once
./build.sh                                 # writes web/cavity.wasm
```

The crate has no dependencies. It exports plain C-ABI functions (`cavity_new`, `cavity_step`, `cavity_u`, ...) that `solver.js` calls directly, so no `wasm-bindgen` or `wasm-pack` is needed.

## Check it against the notebook

```sh
cd web
node validate/check.mjs     # needs python3 with numpy
cargo test --manifest-path cavity/Cargo.toml
```

`check.mjs` runs the notebook's NumPy code (`validate/reference.py`) for 100 and 700 steps on the notebook's 41 × 41 grid. It then compares u, v and p from both the wasm and the JavaScript solvers point by point. The largest difference is about 1e-15, which is round-off. It also prints a rough timing.

Timing on one machine, in ms per step with 50 pressure iterations:

| Grid | NumPy (notebook) | JavaScript | Rust / WASM |
|---|---|---|---|
| 41 × 41 | 1.5 | 0.50 | 0.28 |
| 129 × 129 | 6.7 | 4.6 | 3.0 |

## Files

| File | Purpose |
|---|---|
| `cavity/src/lib.rs` | Rust solver: `build_up_b`, `pressure_poisson` and the time step |
| `solver.js` | `WasmCavity` (loads the module) and `JsCavity`, with the same interface |
| `main.js` | Rendering, tracer particles and controls |
| `index.html` | The page |
| `validate/` | Comparison against the original NumPy code |

## Differences from the notebook

- **Time step:** Δt is the notebook's 0.001 unless the grid or viscosity needs a smaller step to stay stable. On a 129 × 129 grid the notebook's fixed Δt = 0.001 diverges.
- **Grid sizes:** grids up to 257 × 257 are available. The pressure solve is still 50 Jacobi sweeps by default, which converges slowly on fine grids, so raise "Pressure iterations per step" there.
