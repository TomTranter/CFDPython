//! Lid-driven cavity flow, ported from `lessons/14_Step_11.ipynb`.
//!
//! Arrays are stored row-major as `a[j * nx + i]`, matching NumPy's `a[j, i]`:
//! row `j = 0` is the bottom wall (y = 0) and row `j = ny - 1` is the moving lid.
//!
//! The functions at the bottom of this file are a plain C ABI so the module can
//! be loaded with `WebAssembly.instantiate` and no generated JavaScript glue.

pub struct Cavity {
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    pub dt: f64,
    pub rho: f64,
    pub nu: f64,
    pub lid: f64,
    pub nit: usize,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub p: Vec<f64>,
    b: Vec<f64>,
    un: Vec<f64>,
    vn: Vec<f64>,
    pn: Vec<f64>,
    pub steps: u64,
}

impl Cavity {
    /// A cavity of size `lx` x `ly` on an `nx` x `ny` grid, fluid at rest.
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64) -> Self {
        let n = nx * ny;
        Cavity {
            nx,
            ny,
            dx: lx / (nx - 1) as f64,
            dy: ly / (ny - 1) as f64,
            dt: 0.001,
            rho: 1.0,
            nu: 0.1,
            lid: 1.0,
            nit: 50,
            u: vec![0.0; n],
            v: vec![0.0; n],
            p: vec![0.0; n],
            b: vec![0.0; n],
            un: vec![0.0; n],
            vn: vec![0.0; n],
            pn: vec![0.0; n],
            steps: 0,
        }
    }

    pub fn reset(&mut self) {
        for a in [&mut self.u, &mut self.v, &mut self.p, &mut self.b] {
            a.fill(0.0);
        }
        self.steps = 0;
    }

    /// Source term of the pressure Poisson equation (notebook: `build_up_b`).
    fn build_up_b(&mut self) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let (rho, dt) = (self.rho, self.dt);
        let (u, v, b) = (&self.u, &self.v, &mut self.b);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let k = j * nx + i;
                let dudx = (u[k + 1] - u[k - 1]) / (2.0 * dx);
                let dudy = (u[k + nx] - u[k - nx]) / (2.0 * dy);
                let dvdx = (v[k + 1] - v[k - 1]) / (2.0 * dx);
                let dvdy = (v[k + nx] - v[k - nx]) / (2.0 * dy);
                b[k] = rho
                    * (1.0 / dt * (dudx + dvdy) - dudx * dudx - 2.0 * (dudy * dvdx) - dvdy * dvdy);
            }
        }
    }

    /// `nit` Jacobi sweeps of the pressure Poisson equation (notebook: `pressure_poisson`).
    fn pressure_poisson(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let (dx2, dy2) = (self.dx * self.dx, self.dy * self.dy);
        let denom = 2.0 * (dx2 + dy2);
        let (p, pn, b) = (&mut self.p, &mut self.pn, &self.b);
        for _ in 0..self.nit {
            pn.copy_from_slice(p);
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let k = j * nx + i;
                    p[k] = ((pn[k + 1] + pn[k - 1]) * dy2 + (pn[k + nx] + pn[k - nx]) * dx2)
                        / denom
                        - dx2 * dy2 / denom * b[k];
                }
            }
            // Same order as the notebook, which decides the corner values.
            for j in 0..ny {
                p[j * nx + nx - 1] = p[j * nx + nx - 2]; // dp/dx = 0 at x = 2
            }
            for i in 0..nx {
                p[i] = p[nx + i]; // dp/dy = 0 at y = 0
            }
            for j in 0..ny {
                p[j * nx] = p[j * nx + 1]; // dp/dx = 0 at x = 0
            }
            for i in 0..nx {
                p[(ny - 1) * nx + i] = 0.0; // p = 0 at y = 2
            }
        }
    }

    /// Advance one time step (one iteration of the notebook's `cavity_flow` loop).
    pub fn step(&mut self) {
        self.un.copy_from_slice(&self.u);
        self.vn.copy_from_slice(&self.v);

        self.build_up_b();
        self.pressure_poisson();

        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let (dt, rho, nu) = (self.dt, self.rho, self.nu);
        let (u, v, p, un, vn) = (&mut self.u, &mut self.v, &self.p, &self.un, &self.vn);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let k = j * nx + i;
                let (uc, vc) = (un[k], vn[k]);
                u[k] = uc
                    - uc * dt / dx * (uc - un[k - 1])
                    - vc * dt / dy * (uc - un[k - nx])
                    - dt / (2.0 * rho * dx) * (p[k + 1] - p[k - 1])
                    + nu * (dt / (dx * dx) * (un[k + 1] - 2.0 * uc + un[k - 1])
                        + dt / (dy * dy) * (un[k + nx] - 2.0 * uc + un[k - nx]));
                v[k] = vc
                    - uc * dt / dx * (vc - vn[k - 1])
                    - vc * dt / dy * (vc - vn[k - nx])
                    - dt / (2.0 * rho * dy) * (p[k + nx] - p[k - nx])
                    + nu * (dt / (dx * dx) * (vn[k + 1] - 2.0 * vc + vn[k - 1])
                        + dt / (dy * dy) * (vn[k + nx] - 2.0 * vc + vn[k - nx]));
            }
        }

        for i in 0..nx {
            u[i] = 0.0;
            v[i] = 0.0;
            v[(ny - 1) * nx + i] = 0.0;
        }
        for j in 0..ny {
            u[j * nx] = 0.0;
            u[j * nx + nx - 1] = 0.0;
            v[j * nx] = 0.0;
            v[j * nx + nx - 1] = 0.0;
        }
        for i in 0..nx {
            u[(ny - 1) * nx + i] = self.lid; // the moving lid
        }

        self.steps += 1;
    }
}

// ---- C ABI for JavaScript ------------------------------------------------

/// # Safety
/// Every `c` below must be a pointer returned by `cavity_new` and not yet freed.
#[no_mangle]
pub extern "C" fn cavity_new(nx: usize, ny: usize, lx: f64, ly: f64) -> *mut Cavity {
    Box::into_raw(Box::new(Cavity::new(nx, ny, lx, ly)))
}

#[no_mangle]
pub unsafe extern "C" fn cavity_free(c: *mut Cavity) {
    drop(Box::from_raw(c));
}

#[no_mangle]
pub unsafe extern "C" fn cavity_reset(c: *mut Cavity) {
    (*c).reset();
}

#[no_mangle]
pub unsafe extern "C" fn cavity_set_params(c: *mut Cavity, dt: f64, rho: f64, nu: f64, lid: f64, nit: usize) {
    let c = &mut *c;
    c.dt = dt;
    c.rho = rho;
    c.nu = nu;
    c.lid = lid;
    c.nit = nit;
}

#[no_mangle]
pub unsafe extern "C" fn cavity_step(c: *mut Cavity, n: u32) {
    let c = &mut *c;
    for _ in 0..n {
        c.step();
    }
}

#[no_mangle]
pub unsafe extern "C" fn cavity_steps(c: *const Cavity) -> f64 {
    (*c).steps as f64
}

#[no_mangle]
pub unsafe extern "C" fn cavity_u(c: *const Cavity) -> *const f64 {
    (*c).u.as_ptr()
}

#[no_mangle]
pub unsafe extern "C" fn cavity_v(c: *const Cavity) -> *const f64 {
    (*c).v.as_ptr()
}

#[no_mangle]
pub unsafe extern "C" fn cavity_p(c: *const Cavity) -> *const f64 {
    (*c).p.as_ptr()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lid_and_walls_hold_boundary_values() {
        let mut c = Cavity::new(21, 21, 2.0, 2.0);
        for _ in 0..20 {
            c.step();
        }
        let (nx, ny) = (21, 21);
        for i in 1..nx - 1 {
            assert_eq!(c.u[(ny - 1) * nx + i], 1.0);
            assert_eq!(c.u[i], 0.0);
            assert_eq!(c.p[(ny - 1) * nx + i], 0.0);
        }
        assert!(c.u.iter().chain(&c.v).chain(&c.p).all(|x| x.is_finite()));
    }

    #[test]
    fn reset_returns_to_rest() {
        let mut c = Cavity::new(11, 11, 2.0, 2.0);
        c.step();
        c.reset();
        assert!(c.u.iter().chain(&c.v).chain(&c.p).all(|&x| x == 0.0));
        assert_eq!(c.steps, 0);
    }
}
