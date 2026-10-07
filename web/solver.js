// Two interchangeable cavity-flow solvers with the same interface:
//   setParams({dt, rho, nu, lid, nit}), step(n), reset(), fields() -> {u, v, p}, free()
// WasmCavity calls the Rust module in cavity/; JsCavity is the same scheme in
// plain JavaScript, kept for speed comparison and as a fallback.

export async function loadWasm(bytes) {
  const { instance } = await WebAssembly.instantiate(bytes, {});
  return instance.exports;
}

export class WasmCavity {
  constructor(wasm, nx, ny, lx = 2, ly = 2) {
    this.w = wasm;
    this.nx = nx;
    this.ny = ny;
    this.ptr = wasm.cavity_new(nx, ny, lx, ly);
  }
  setParams({ dt, rho, nu, lid, nit }) {
    this.w.cavity_set_params(this.ptr, dt, rho, nu, lid, nit);
  }
  step(n = 1) {
    this.w.cavity_step(this.ptr, n);
  }
  reset() {
    this.w.cavity_reset(this.ptr);
  }
  get steps() {
    return this.w.cavity_steps(this.ptr);
  }
  // Views into wasm memory. Rebuild them after any call that could grow memory.
  fields() {
    const buf = this.w.memory.buffer;
    const n = this.nx * this.ny;
    return {
      u: new Float64Array(buf, this.w.cavity_u(this.ptr), n),
      v: new Float64Array(buf, this.w.cavity_v(this.ptr), n),
      p: new Float64Array(buf, this.w.cavity_p(this.ptr), n),
    };
  }
  free() {
    this.w.cavity_free(this.ptr);
  }
}

export class JsCavity {
  constructor(nx, ny, lx = 2, ly = 2) {
    this.nx = nx;
    this.ny = ny;
    this.dx = lx / (nx - 1);
    this.dy = ly / (ny - 1);
    const n = nx * ny;
    for (const k of ["u", "v", "p", "b", "un", "vn", "pn"]) this[k] = new Float64Array(n);
    this.setParams({ dt: 0.001, rho: 1, nu: 0.1, lid: 1, nit: 50 });
    this.steps = 0;
  }
  setParams({ dt, rho, nu, lid, nit }) {
    Object.assign(this, { dt, rho, nu, lid, nit });
  }
  reset() {
    for (const k of ["u", "v", "p", "b"]) this[k].fill(0);
    this.steps = 0;
  }
  fields() {
    return { u: this.u, v: this.v, p: this.p };
  }
  free() {}

  step(count = 1) {
    const { nx, ny, dx, dy, dt, rho, nu, nit, u, v, p, b, un, vn, pn } = this;
    const dx2 = dx * dx, dy2 = dy * dy, denom = 2 * (dx2 + dy2);
    for (let s = 0; s < count; s++) {
      un.set(u);
      vn.set(v);

      for (let j = 1; j < ny - 1; j++) {
        for (let i = 1; i < nx - 1; i++) {
          const k = j * nx + i;
          const dudx = (u[k + 1] - u[k - 1]) / (2 * dx);
          const dudy = (u[k + nx] - u[k - nx]) / (2 * dy);
          const dvdx = (v[k + 1] - v[k - 1]) / (2 * dx);
          const dvdy = (v[k + nx] - v[k - nx]) / (2 * dy);
          b[k] = rho * (1 / dt * (dudx + dvdy) - dudx * dudx - 2 * (dudy * dvdx) - dvdy * dvdy);
        }
      }

      for (let q = 0; q < nit; q++) {
        pn.set(p);
        for (let j = 1; j < ny - 1; j++) {
          for (let i = 1; i < nx - 1; i++) {
            const k = j * nx + i;
            p[k] = ((pn[k + 1] + pn[k - 1]) * dy2 + (pn[k + nx] + pn[k - nx]) * dx2) / denom
              - dx2 * dy2 / denom * b[k];
          }
        }
        for (let j = 0; j < ny; j++) p[j * nx + nx - 1] = p[j * nx + nx - 2];
        for (let i = 0; i < nx; i++) p[i] = p[nx + i];
        for (let j = 0; j < ny; j++) p[j * nx] = p[j * nx + 1];
        for (let i = 0; i < nx; i++) p[(ny - 1) * nx + i] = 0;
      }

      for (let j = 1; j < ny - 1; j++) {
        for (let i = 1; i < nx - 1; i++) {
          const k = j * nx + i;
          const uc = un[k], vc = vn[k];
          u[k] = uc
            - uc * dt / dx * (uc - un[k - 1])
            - vc * dt / dy * (uc - un[k - nx])
            - dt / (2 * rho * dx) * (p[k + 1] - p[k - 1])
            + nu * (dt / dx2 * (un[k + 1] - 2 * uc + un[k - 1])
              + dt / dy2 * (un[k + nx] - 2 * uc + un[k - nx]));
          v[k] = vc
            - uc * dt / dx * (vc - vn[k - 1])
            - vc * dt / dy * (vc - vn[k - nx])
            - dt / (2 * rho * dy) * (p[k + nx] - p[k - nx])
            + nu * (dt / dx2 * (vn[k + 1] - 2 * vc + vn[k - 1])
              + dt / dy2 * (vn[k + nx] - 2 * vc + vn[k - nx]));
        }
      }

      for (let i = 0; i < nx; i++) {
        u[i] = 0;
        v[i] = 0;
        v[(ny - 1) * nx + i] = 0;
      }
      for (let j = 0; j < ny; j++) {
        u[j * nx] = 0;
        u[j * nx + nx - 1] = 0;
        v[j * nx] = 0;
        v[j * nx + nx - 1] = 0;
      }
      for (let i = 0; i < nx; i++) u[(ny - 1) * nx + i] = this.lid;
      this.steps++;
    }
  }
}
