# Handoff notes: browser CFD exploration

Notes from a Claude Code session on 2026-10-07, written so the work can be picked up in a new session.
Branch: `claude/amazing-mccarthy-n3pu0m` on `TomTranter/CFDPython`.

## What happened, in order

1. **Synced with upstream.** This fork was 104 commits behind `barbagroup/CFDPython` `master`. The branch was fast-forwarded to upstream (`92c8220`); no local changes were lost. The fork's own `master` was not touched.
2. **Discussed porting to a faster language for the browser.** Recommendation: Rust compiled to WebAssembly, with a light JavaScript front end. The notebooks use tiny grids (41 × 41, at most 700 steps), so the gain is interactivity and larger grids rather than raw speed. WebGPU is the next step if the pressure solve on large grids becomes the bottleneck.
3. **Built a prototype of Step 11 (lid-driven cavity)** in this folder. See [`README.md`](README.md) for how to run, rebuild and validate it.
4. **Found a crash at low viscosity.** Not yet fixed; details below.
5. **Researched existing tools** instead of reinventing the wheel: OpenFOAM in the browser, OpenFOAM-based cloud services, desktop interfaces for OpenFOAM, and HELYX. Summary below.
6. **Idea: use this as a showcase for Modlit (modlit.io).** Not started. The session couldn't read modlit.io because the cloud environment's network policy blocks that domain, and a web search found nothing about it.

## The prototype

| Path | What it is |
|---|---|
| `cavity/src/lib.rs` | Rust port of the notebook's `build_up_b`, `pressure_poisson` and time step. No dependencies; plain C-ABI exports, so no wasm-bindgen. |
| `cavity.wasm` | The compiled module (27 KB), committed so the page works without Rust installed. |
| `solver.js` | `WasmCavity` and `JsCavity`, two solvers with the same interface. |
| `main.js`, `index.html` | The interactive page: pressure, speed or vorticity colour map, tracer particles, arrows, sliders, and an engine switch. |
| `validate/` | `reference.py` holds the notebook's NumPy code; `check.mjs` compares both solvers against it. |

- **Accuracy.** Both solvers match the notebook to about 1e-15 (round-off) after 100 and 700 steps on the 41 × 41 grid. `cargo test` passes (two boundary-condition and reset tests).
- **Speed** (ms per step, 50 pressure iterations, Node and this container):

  | Grid | NumPy | JavaScript | Rust / WASM |
  |---|---|---|---|
  | 41 × 41 | 1.5 | 0.50 | 0.28 |
  | 129 × 129 | 6.7 | 4.6 | 3.0 |

- **Time step.** The page shrinks Δt below the notebook's 0.001 when the grid or viscosity needs it: `dt = min(0.001, 0.2 h²/ν, 0.25 h/U)`. The notebook's own code diverges at 129 × 129 with Δt = 0.001.
- **Live copy.** Published as a private Claude artifact: https://claude.ai/artifact/AwxNJmRpmyhpSwbP4f5fE9. The artifact copy is built from `index.html` with the `<html>`, `<head>` and `<body>` wrapper removed. Nobody checked whether WebAssembly loads in the artifact viewer; if it doesn't, the page falls back to the JavaScript solver and says so.

## Open bug: crash at low viscosity

The user reported that the simulation crashes (the field fills with NaN or garbage) when the viscosity slider is low. A reproduction script was about to run when the session moved on, so **the cause is not confirmed.** Likely causes, in order:

1. **Downwind differencing.** The notebook uses backward differences, such as `u * (u[i] - u[i-1])`, whatever the sign of the velocity. Where the flow runs the other way (the return flow along the bottom and sides of the cavity), that becomes a downwind scheme, which is unstable unless diffusion damps it. At low ν, the grid Reynolds number `U·h/ν` is large (10 at 41 × 41 with ν = 0.005), so nothing damps it.
2. **Pressure under-resolved.** 50 Jacobi sweeps per step leave the pressure far from converged, especially on fine grids, so the velocity field isn't kept close to divergence-free.
3. **Advective time-step limit uses the lid speed U.** Internal velocities rarely exceed U, so this is the least likely cause.

Suggested fix:
- Add a sign-aware upwind option in both `lib.rs` and `solver.js`. Keep the notebook's scheme available, so `validate/check.mjs` can still compare against NumPy exactly.
- Make upwind the page's default.
- Detect non-finite values after each frame, pause, and tell the user to raise viscosity or reset.
- Consider limiting the slider's lower end by grid size.

Reproduce with a loop like this before fixing:
```js
// For each grid n in [41, 81, 129] and nu in [0.05, 0.02, 0.01, 0.005]:
// run WasmCavity with the page's dt rule and report the step at which max |u| exceeds 10 or goes non-finite.
```

## Research summary: OpenFOAM and the browser

Mostly from the model's background knowledge, with light web searching. Check the claims marked "(unverified)" before relying on them.

- **No major open-source CFD package runs its solver in the browser.** That covers OpenFOAM, SU2, FEniCS and similar. The obstacles are huge C++ codebases, solvers loaded as shared libraries at runtime, file-based case setup, MPI, and memory limits. A single-core Emscripten build of OpenFOAM is conceivable but would be a project of its own.
- **What does run in the browser:**
  - WebGPU and WebGL lattice-Boltzmann solvers, for example Polymère.
  - Small hand-written solvers like this one.
  - Post-processing with vtk.js or VTK's WebAssembly build.
  - The original notebooks themselves, through JupyterLite or Pyodide (unverified; the Numba lesson won't work, because Numba isn't available in Pyodide).
- **OpenFOAM-based cloud services:**
  - SimScale and AirShaper: full SaaS, with OpenFOAM hidden from the user.
  - CFD Direct's CFDDFC: an AWS or Azure machine image with OpenFOAM v14, ParaView and FreeCAD; "Web CFDDFC" adds a remote desktop in the browser.
  - Rescale, Inductiva and CloudHPC: on-demand computing.
- **Web CFDDFC isn't open source** (unverified; cfd.direct is blocked from this environment). It bundles open-source software: OpenFOAM is GPL v3, ParaView is BSD and FreeCAD is LGPL. The image and its web desktop are CFD Direct's product. An equivalent can be built from open source: a cloud Ubuntu machine with OpenFOAM and ParaView, plus Apache Guacamole or noVNC.
- **Desktop interfaces for OpenFOAM:**
  - FreeCAD with the CfdOF add-on: free, and the best free option.
  - HELYX-OS: open source, but likely unmaintained.
  - SimFlow: free tier.
  - ENGYS HELYX, ESI Visual-CFD, CFD Support TCFD and CastNet: commercial.
  - CFDTool / FEATool: runs in MATLAB.
- **HELYX** is ENGYS's commercial OpenFOAM-based package. It has an extended solver, its own mesher derived from snappyHexMesh, and add-on modules for optimisation and multiphase flow. It is aimed at automotive, building ventilation, marine and turbomachinery work.

**Conclusion from the discussion.** The user does not want to compete in the CFD-tool market. The interest is in this work as a showcase.

## Next steps for the new session

1. **Learn what Modlit is.** Read modlit.io (allow the domain in the environment's network settings) or get a description from the user. Then decide what the showcase should demonstrate. Options discussed:
   - Notebook to live app.
   - Performance and WebAssembly.
   - Wrapping existing solvers.
2. **Fix the low-viscosity crash** (see above) before showing the demo to anyone.
3. **Possibly extend the demo across lessons.** The 12 lessons build on each other, so the showcase could be a series rather than a single demo. The licences allow reuse: code is BSD-3, text is CC-BY 4.0.
4. **Possibly add a JupyterLite build** of the lessons, if running the notebooks in the browser unchanged is enough for the showcase.
