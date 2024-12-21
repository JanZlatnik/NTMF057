/*

RADIALSE MODULE

Contains: Procedures to compute eigenstates of radial Schrödinger equation

Last revision: 05/12/2024 by JZ

*/
use nalgebra::{DMatrix, DVector};
use crate::core::console;
use crate::math::{runge_kutta, Grid, Integrable, ODESetup, RKMethod, ODE};
use crate::math::INTMethod::Simpson;


// Shooting method \\
pub fn calculate_wavefunction <F> (
    grid: &Grid,
    potential: F,
    l: usize,
    mass: f64,
    energy: f64
) -> Vec<f64>
where
    F: Fn(f64) -> f64 + Clone + 'static
{
    let ini = vec![0.0,f64::EPSILON];
    let schrodinger = ODE::new(2,move |state: &[f64],r| {
        if state[0] == 0.0 {
            return 0.0;
        }
        ((l*(l+1)) as f64 / r.powi(2) + 2.0 * mass * potential(r) - 2.0 * mass * energy) * state[0]
    });
    let psi = runge_kutta(
        &*schrodinger,
        &ini,
        &grid,
        RKMethod::RK4
    );
    let wavefunction = psi[0].clone();
    wavefunction
}

fn asymptotic(psi: Vec<f64>) -> f64 {
    let x = psi[psi.len() - 1];
    x.tanh()
}

pub fn shooting_method <F>(
    grid: &Grid,
    potential: F,
    l: usize,
    mass: f64,
    ne: usize
) -> (DVector<f64>, DMatrix<f64>, DVector<f64>, DVector<f64>)
where
    F: Fn(f64) -> f64 + Clone + 'static,{

    let effective_potential: Vec<f64> = grid.points
        .iter()
        .map(|&r| potential(r) + (l * (l + 1)) as f64 / (2.0 * mass * r.powi(2)))
        .collect();

    let emin = effective_potential
        .iter()
        .filter(|&&x| x.is_finite())
        .copied()
        .reduce(f64::min)
        .unwrap_or(-f64::INFINITY);

    let emax = 0.0; //potential(10.0 * grid.end) + (l * (l + 1)) as f64 / (2.0 * mass * (10.0 * grid.end).powi(2));
    console(&format!("Searching for bound states in energy range [{emin},{emax}] au"));

    let energies = Grid::new_n(emin,emax,ne);
    let sign_function: Vec<f64> = energies.points.iter().map(|&e| { asymptotic(calculate_wavefunction(&grid, potential.clone(), l, mass, e)) }).collect();

    let mut bound_energies: Vec<f64> = Vec::new();

    for (i, window) in sign_function.windows(2).enumerate() {
        if window[0] * window[1] <= 0.0 {
            let mut e_left = energies.points[i];
            let mut e_right = energies.points[i + 1];
            let mut converged = false;

            let tolerance = 1e-14;
            for _ in 0..100 {
                let e_mid = (e_left + e_right) / 2.0;
                let sign_mid = asymptotic(calculate_wavefunction(&grid, potential.clone(), l, mass, e_mid));

                if (e_left - e_right).abs() <= tolerance {
                    bound_energies.push(e_mid);
                    converged = true;
                    break;
                }

                if sign_mid * window[0] <= 0.0 {
                    e_right = e_mid;
                } else {
                    e_left = e_mid;
                }
            }

            if !converged {
                console(&format!("[ERROR]: Bisection did not converge at energy {} au", e_left));
            }
        }
    }
    let final_energies = DVector::from_vec(bound_energies);
    let mut eigenfunctions = DMatrix::zeros(final_energies.len(),grid.n);

    for (i,&e) in final_energies.iter().enumerate() {
        let mut psi = DVector::from_vec(calculate_wavefunction(&grid, potential.clone(), l, mass, e));
        let last_classical_turning_point = effective_potential.iter().rposition(|&pot| pot <= e).unwrap_or(0);
        let mut cutoff_index = grid.n - 1;
        for j in (last_classical_turning_point..grid.n-1) {
            if psi[j+1].abs() > psi[j].abs() {
                cutoff_index = j;
                break;
            }
        }
        for j in cutoff_index+1..grid.n {
            psi[j] = 0.0;
        }
        let norm = psi.square_integrate(grid.d, Simpson);
        let normalized_psi =  psi * (1.0 / norm.sqrt());

        eigenfunctions.row_mut(i).copy_from_slice(&normalized_psi.as_slice());
    }
    (final_energies,eigenfunctions,DVector::from_vec(energies.points.clone()),DVector::from_vec(sign_function))
}

