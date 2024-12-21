/*

MATH MODULE

Contains: Basic mathematical functions and operations

Last revision: 21/12/2024 by JZ

*/
use std::fmt::Debug;
use nalgebra::{Complex, DMatrix, DVector};
use num::traits::{Num, Zero};
use std::ops::{Add, Div, Mul};
use num_traits::{Float, NumCast, Signed};
use num::complex::ComplexFloat;
use crate::core::{HBAR, PI};


// Grid definition \\
#[derive(Clone)]
pub struct Grid {
    pub points: Vec<f64>,
    pub d:      f64,
    pub n:      usize,
    pub start:  f64,
    pub end:    f64,
}

impl Grid {
    pub(crate) fn new_n(start: f64, end: f64, n: usize) -> Grid {
        let d = (end - start)/((n - 1) as f64);
        let points = (0.. n)
            .map(|i| start + i as f64 * d)
            .collect();

        Grid{
            points,
            d,
            n,
            start,
            end
        }
    }

    pub(crate) fn new_d(start: f64, end: f64, spacing: f64) -> Grid {
        let d = spacing;
        let n = ((end - start) / d).ceil() as usize + 1;

        let points = (0..n)
            .map(|i| start + i as f64 * d)
            .collect();


        Grid {
            points,
            d,
            n,
            start,
            end,
        }
    }

    pub(crate) fn tail(&mut self, t:f64) {
        if t == 1.0 {
            return;
        }
        else if t < 0.0 {
            return;
        }
        else if t >  1.0 {
            let new_n = ((self.n as f64) * t).ceil() as usize;
            for i in 0..(new_n-self.n) {
                self.points.push(self.points[self.n-1+i] + self.d);
            }
            self.end = self.points[new_n-1];
            self.n = self.points.len();
        }
        else if t < 1.0 {
            let new_n = ((self.n as f64) * t).ceil() as usize;
            self.points.truncate(new_n);
            self.end = self.points[new_n-1];
            self.n = self.points.len();
        }
    }
}

// Transpose \\
pub trait Transposable<T> {
    fn transpose(&self) -> Vec<Vec<T>>;
}
impl<T: Clone + Zero> Transposable<T> for Vec<Vec<T>> {
    fn transpose(&self) -> Vec<Vec<T>> {
        if self.is_empty() {
            return vec![];
        }

        let row_count = self.len();
        let col_count = self[0].len();

        let mut transposed = vec![vec![T::zero(); row_count]; col_count];

        for (i, row) in self.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                transposed[j][i] = value.clone();
            }
        }

        transposed
    }
}


// Integration \\
#[derive(Clone, Copy)]
pub enum INTMethod {
    Trapezoid,
    Simpson,
    INT4
}
#[derive(Clone, Copy)]
pub enum PrimitiveMethod {
    P2,
    P4
}
pub trait Integrable<T> {
    fn integrate(&self, dx: f64, method: INTMethod) -> T;
    fn square_integrate(&self, dx: f64, method: INTMethod) -> f64;
    fn primitive(&self, dx: f64, method: PrimitiveMethod) -> DVector<T>;
}

trait ComplexOrReal {
    fn abs(&self) -> f64;
}
impl ComplexOrReal for f64 {
    fn abs(&self) -> f64 { num_traits::Float::abs(*self)  }
}
impl<F: Float> ComplexOrReal for Complex<F> {
    fn abs(&self) -> f64 { self.norm().to_f64().unwrap_or(0.0) }
}

impl<T: Num + Copy + Add<Output=T> + Mul<f64, Output=T> + Zero + Default + Debug + 'static> Integrable<T> for DVector<T>
where T: ComplexOrReal {
    fn integrate(&self, dx: f64, method: INTMethod) -> T {
        let n = self.len();
        let mut sum = T::zero();
        match method {
            INTMethod::Trapezoid => {
                if n < 2 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                for i in 1..n-1 {
                    sum = sum + self[i];
                }
                sum = sum + (self[0]+self[n-1]) * 0.5;
            }
            INTMethod::Simpson => {
                if n < 3 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                let m = n - 1 + (n % 2);
                sum = sum + self[0] + self[m-1];
                for i in 1..=(m - 1) / 2 {
                    sum = sum + self[2*i-1] * 4.0
                }
                for i in 1..=(m-1)/2 - 1 {
                    sum = sum + self[2*i] * 2.0
                }
                sum = sum * (1.0 / 3.0);
                if m != n {
                    sum = sum + self[n-1] * (3.0/8.0) + self[n-2] * (19.0/24.0) - self[n-3] * (5.0/24.0) + self[n-4] * (1.0/ 24.0);
                }
            }
            INTMethod::INT4 => {
                if n < 6 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                for i in 3..n-3 {
                    sum = sum + self[i]
                }
                sum = sum + (self[0]+self[n-1]) * (3.0/8.0) + (self[1]+self[n-2]) * (7.0 / 6.0) + (self[2]+self[n-3]) * (23.0 / 24.0)
            }
        }

        sum * dx
    }

    fn square_integrate(&self, dx: f64, method: INTMethod) -> f64 {
        let squared_vector = DVector::from_vec(
            self.iter().map(|&x| x.abs() * x.abs()).collect()
        );
        squared_vector.integrate(dx, method)
    }


    fn primitive(&self, dx: f64, method: PrimitiveMethod) -> DVector<T> {
        let n = self.len();
        let mut primitive_vec =  DVector::from_vec(vec![T::zero(); n]);
        match method {
            PrimitiveMethod::P2 => {
                if n < 2 {
                    panic!("[Primitive Error]: Not enough elements")
                }
                for i in 1..n {
                    primitive_vec[i] = primitive_vec[i-1] + (self[i-1]+self[i])* 0.5  * dx;
            }
            }
            PrimitiveMethod::P4 => {
                if n < 4 {
                    panic!("[Primitive Error]: Not enough elements")
                }
                primitive_vec[1] = (self[0] * (3.0/8.0)  + self[1] * (19.0/24.0) - self[2] * (5.0/24.0) + self[3] * (1.0/24.0))  * dx;
                for i in 2..n-1 {
                    primitive_vec[i] = primitive_vec[i-1] + ((self[i-1]+self[i]) * (13.0/24.0) - (self[i-2]+self[i+1]) * (1.0/24.0)) * dx;
                }
                primitive_vec[n-1] = primitive_vec[n-2] + (self[n-1] * (3.0/8.0) + self[n-2] * (19.0/24.0) - self[n-3] * (5.0/24.0) + self[n-4] * (1.0/24.0)) * dx;
            }

        }
        primitive_vec
    }
}


// Eigenstates of radial Schrödinger equation \\
pub struct Eigenstates {
    pub e: DVector<f64>,
    pub four: DMatrix<f64>,
    pub func: DMatrix<f64>,
}

impl Eigenstates {
    pub(crate) fn fourier_dvr(
        grid: &Grid,
        potential: &dyn Fn(f64) -> f64,
        mass: f64,
        l: usize,
        nfdvr: usize,
        n:  usize,
    ) -> Result<Self,String> {
        if n > nfdvr {
            return Err("[FDVR Error]: n is greater than nfdvr".into());
        }

        let effective_potential = |r: f64| -> f64  {potential(r) + (l*(l+1)) as f64 / (2.0 * mass * r.powi(2))};
        let length = (grid.end-grid.start).abs();
        let mut x_optor = DMatrix::<f64>::zeros(nfdvr, nfdvr);

        for i in 2..=nfdvr{
            for j in ((i % 2) + 1..i).step_by(2) {
                let value = - 8.0 * length * (i as f64) * (j as f64) / (PI.powi(2) * (i.pow(2) as f64 - j.pow(2) as f64).powi(2));
                x_optor[(i-1,j-1)] = value;
                x_optor[(j-1,i-1)] = value;
            }
        }

        let x_eig_res = x_optor.clone().symmetric_eigen();
        let mut x_eig = x_eig_res.eigenvalues.add_scalar((grid.start+grid.end)* 0.5);
        let x_eig_vectors = x_eig_res.eigenvectors;

        let mut h_matrix = DMatrix::<f64>::zeros(nfdvr, nfdvr);
        for i in 0..nfdvr{
            for j in 0..=i{
                let mut h_value = 0.0;
                for k in 0..nfdvr{
                    h_value += x_eig_vectors[(i,k)] * effective_potential(x_eig[k]) * x_eig_vectors[(j, k)];
                }
                if i == j {
                    h_value += ((i+1).pow(2) as f64) * PI.powi(2) * HBAR.powi(2) / (2.0 * mass * length.powi(2));
                }
                h_matrix[(i,j)] = h_value;
                h_matrix[(j,i)] = h_value;
            }
        }

        let eig_res = h_matrix.symmetric_eigen();
        let eigenvalues = eig_res.eigenvalues.as_slice();
        let eigenvectors = eig_res.eigenvectors.clone();
        let mut indices: Vec<usize> = (0..eigenvalues.len()).collect();
        indices.sort_by(|&i, &j| eigenvalues[i].partial_cmp(&eigenvalues[j]).unwrap());
        let sorted_eigenvalues: Vec<f64> = indices.iter().map(|&i| eigenvalues[i]).collect();
        let sorted_eigenvectors: DMatrix<f64> = DMatrix::from_fn(eigenvectors.nrows(), eigenvectors.ncols(), |row, col| {
            eigenvectors[(row, indices[col])]
        });
        let e = DVector::from_vec(sorted_eigenvalues[0..n].to_vec());
        let four: DMatrix<f64> = sorted_eigenvectors.columns(0, n).into();

        let mut func = DMatrix::<f64>::zeros(n, grid.n);
        for k in 0..n {
            for i in 0..nfdvr {
                for j in 0..grid.n {
                    func[(k,j)] += four[(i,k)] * (2.0 / length).sqrt()
                        * (PI * (grid.points[j]-grid.start)*(i+1) as f64 / length).sin();
                }
            }
        }

        Ok(Self{
            e,
            four,
            func
        })

    }
}



// Riccati-Bessel and Riccati-Neumann function definitions for real arguments \\
pub fn riccati_bessel(l: usize, x: f64) -> f64 {
    match l {
        0 => x.sin(),
        1 => x.sin() / x - x.cos(),
        _ => {
            let mut j0 = x.sin();
            let mut j1 = x.sin() / x - x.cos();
            let mut jl = 0.0;
            for i in 2..=l {
                jl = ((2 * i - 1) as f64 / x) * j1 - j0;
                j0 = j1;
                j1 = jl;
            }
            jl
        }
    }
}

pub fn riccati_neumann(l: usize, x: f64) -> f64 {
    match l {
        0 => x.cos(),
        1 => x.cos() / x + x.sin(),
        _ => {
            let mut n0 = x.cos();
            let mut n1 = x.cos() / x + x.sin();
            let mut nl = 0.0;
            for i in 2..=l {
                nl = ((2 * i - 1) as f64 / x) * n1 - n0;
                n0 = n1;
                n1 = nl;
            }
            nl
        }
    }
}



// Runge - Kutta \\
#[derive(Clone, Copy)]
pub enum RKMethod {
    RK1,
    RK2,
    RK4
}

pub trait GeneralODE<T> {
    fn order(&self) -> usize;
    fn right_side(&self, state: &[T], x: f64) -> T;
}

pub trait ODESetup<T> {
    fn new<F>(order: usize, right_side_fn: F) -> Box<dyn GeneralODE<T>>
    where
        F: Fn(&[T],f64) -> T + 'static;
}

pub struct ODE;

impl<T> ODESetup<T> for ODE
where
    T:  NumCast + 'static,
{
    fn new<F>(order: usize, right_side_fn: F) -> Box<dyn GeneralODE<T>>
    where
        F: Fn(&[T],f64) -> T + 'static
    {
        struct DynamicODE<F> {
            order: usize,
            right_side_fn: F,
        }
        impl<F,T> GeneralODE<T> for DynamicODE<F>
        where
            F: Fn(&[T],f64) -> T,
            T: NumCast,
        {
            fn order(&self) -> usize {
                self.order
            }
            fn right_side(&self, state: &[T], x: f64) -> T {
                (self.right_side_fn)(state,x)
            }
        }
        Box::new(DynamicODE {
            order,
            right_side_fn
        })
    }
}

pub fn runge_kutta<T> (
    ode: &dyn GeneralODE<T>,
    initial_conditions: &[T],
    grid: &Grid,
    method: RKMethod,
) -> Vec<Vec<T>>
where
    T: NumCast + Add<Output = T> + Mul<Output = T> + Div<Output = T> + Clone + Zero + Copy
{
    let n = ode.order();
    let mut solution = Vec::with_capacity(grid.n);
    let mut current_state = initial_conditions.to_vec();
    solution.push(current_state.clone());

    match method {
        RKMethod::RK1 => {
            for i in 0..grid.n-1 {
                let x: f64 = grid.points[i];
                let h = T::from(grid.d).unwrap();
                let mut next_state = current_state.clone();

                let mut k1 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k1[j] = current_state[j+1];
                }
                k1[n-1] = ode.right_side(&current_state, x);


                for j in 0..n {
                    next_state[j] = next_state[j] + h * k1[j];
                }

                current_state = next_state;
                solution.push(current_state.clone());
            }
        }

        RKMethod::RK2 => {
            for i in 0..grid.n-1 {
                let x: f64 = grid.points[i];
                let h = T::from(grid.d).unwrap();
                let mut next_state = current_state.clone();

                let mut k1 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k1[j] = current_state[j+1];
                }
                k1[n-1] = ode.right_side(&current_state, x);

                let mut intermediate_state = current_state.clone();
                for j in 0..n {
                    intermediate_state[j] = intermediate_state[j] + h * k1[j] / T::from(2.0).unwrap();
                }

                let mut k2 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k2[j] = intermediate_state[j+1];
                }
                k2[n-1] = ode.right_side(&intermediate_state, x + grid.d / 2.0);

                for j in 0..n {
                    next_state[j] = next_state[j] + h * k2[j];
                }

                current_state = next_state;
                solution.push(current_state.clone());
            }
        }

        RKMethod::RK4 => {
            for i in 0..grid.n-1 {
                let x: f64 = grid.points[i];
                let h = T::from(grid.d).unwrap();
                let mut next_state = current_state.clone();

                let mut k1 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k1[j] = current_state[j+1];
                }
                k1[n-1] = ode.right_side(&current_state, x);

                let mut intermediate_state1 = current_state.clone();
                for j in 0..n {
                    intermediate_state1[j] = intermediate_state1[j] + h * k1[j] / T::from(2.0).unwrap();
                }
                let mut k2 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k2[j] = intermediate_state1[j+1];
                }
                k2[n-1] = ode.right_side(&intermediate_state1, x + grid.d / 2.0);

                let mut intermediate_state2 = current_state.clone();
                for j in 0..n {
                    intermediate_state2[j] = intermediate_state2[j] + h * k2[j] / T::from(2.0).unwrap();
                }
                let mut k3 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k3[j] = intermediate_state2[j+1];
                }
                k3[n-1] = ode.right_side(&intermediate_state2, x + grid.d / 2.0);

                let mut intermediate_state3 = current_state.clone();
                for j in 0..n {
                    intermediate_state3[j] = intermediate_state3[j] + h * k3[j];
                }
                let mut k4 = vec![T::zero(); n];
                for j in 0..n-1 {
                    k4[j] = intermediate_state3[j+1];
                }
                k4[n-1] = ode.right_side(&intermediate_state3, x + grid.d);

                for j in 0..n {
                    next_state[j] = next_state[j] + (k1[j] + T::from(2.0).unwrap() * k2[j] + T::from(2.0).unwrap() * k3[j] + k4[j]) * h / T::from(6.0).unwrap();
                }

                current_state = next_state;
                solution.push(current_state.clone());
            }
        }

    }
    solution.transpose()
}



// Natural cubic spline \\
#[derive(Clone)]
struct NaturalCubicSpline {
    x_values: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl NaturalCubicSpline {
    fn new(x_val: Vec<f64>, y_val: Vec<f64>) -> Self {
        assert_eq!(x_val.len(), y_val.len(), "[ERROR] x and y values do NOT have the same length.");
        let n = x_val.len() - 1;

        let a = y_val.clone();
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];

        let h: Vec<f64> = x_val.windows(2).map(|w| w[1] - w[0]).collect();

        let mut alpha = vec![0.0; n - 1];
        for i in 1..n {
            alpha[i - 1] = (a[i + 1] - a[i]) * 3.0 / h[i] - (a[i] - a[i - 1]) * 3.0 / h[i - 1];
        }

        let mut c = vec![0.0; n + 1];
        let mut l = vec![0.0; n + 1];
        let mut mu = vec![0.0; n + 1];
        let mut z = vec![0.0; n + 1];

        l[0] = 1.0;
        for i in 1..n {
            l[i] = 2.0 * (x_val[i + 1] - x_val[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n] = 1.0;
        for i in (0..n).rev() {
            c[i] = z[i] - mu[i] * c[i + 1];
            b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }

        NaturalCubicSpline {
            x_values: x_val,
            a,
            b,
            c,
            d,
        }
    }

    fn evaluate(&self, x: f64) -> f64 {
        let n = self.x_values.len() - 1;

        if x < self.x_values[0] || x > self.x_values[n] {
            return f64::NAN;
        }

        let j = match self.x_values.windows(2)
            .position(|w| x >= w[0] && x <= w[1]) {
            Some(index) => index,
            None => return f64::NAN,
        };

        let dx = x - self.x_values[j];
        self.a[j] + self.b[j] * dx + self.c[j] * dx.powi(2) + self.d[j] * dx.powi(3)
    }
}

pub fn create_natural_cubic_spline(x_val: Vec<f64>, y_val: Vec<f64>) -> impl Fn(f64) -> f64 + Clone {
    let spline = NaturalCubicSpline::new(x_val, y_val);
    move |x| spline.evaluate(x)
}
