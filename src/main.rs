use std::fs;
use std::fs::File;
use crate::core::{console, Settings, PHYS_H0, PHYS_ME, PHYS_U};
use crate::math::{create_natural_cubic_spline, runge_kutta, Eigenstates, Grid, Integrable, ODESetup, ODE};
use std::io::{Write};
use nalgebra::{ComplexField};
use crate::math::RKMethod::{RK1, RK2, RK4};
use crate::radialSE::{shooting_method};

mod math;
mod radialSE;
mod core;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    fs::create_dir_all("Task 1")?;
    fs::create_dir_all("Task 3")?;
    fs::create_dir_all("Test")?;

    let settings = Settings::initialize("settings.toml")?;
    let r = Grid::new_n(settings.rmin,settings.rmax,settings.nr);
    let m = 14.007 * PHYS_U / PHYS_ME / 2.0;

    let n2 = |r: f64| -> f64 {
        let v0 = 0.75102;
        let alpha = 1.15350;
        let r0 = 2.01943;
        v0 * ((-2.0*alpha*(r-r0)).exp() - 2.0*(-alpha*(r-r0)).exp())
    };

    // Potential output \\
    let mut file = File::create("Test/N2_potential.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [eV]")?;
    for i in 0..r.n {
        write!(file, "{:20.12E}", r.points[i])?;
        write!(file, "{:20.12E}", n2(r.points[i])*PHYS_H0)?;
        writeln!(file)?;
    }
    console("Data successfully written to Test/N2_potential.txt");

    // Runge Kutta test on LHO \\
    let mut file = File::create("Task 1/RK_test.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}{:>20}{:>20}", "#", "n", "RK1", "RK2", "RK4")?;
    let ini = vec![1.0,0.0];
    let lho_ode = ODE::new(2, |state: &[f64],x| {-state[0]});
    for i in 1..17 {
        let mut n = 50;
        for _ in 0..i { n *= 2 }
        let x = Grid::new_n(0.0,20.0,n);
        let sol_1 = runge_kutta(&*lho_ode, &ini, &x, RK1);
        let sol_2 = runge_kutta(&*lho_ode, &ini, &x, RK2);
        let sol_4 = runge_kutta(&*lho_ode, &ini, &x, RK4);
        write!(file, "{:20.0}", n)?;
        write!(file, "{:20.12E}", (sol_1[0][x.n-1] - x.points[x.n-1].cos()).abs())?;
        write!(file, "{:20.12E}", (sol_2[0][x.n-1] - x.points[x.n-1].cos()).abs())?;
        write!(file, "{:20.12E}", (sol_4[0][x.n-1] - x.points[x.n-1].cos()).abs())?;
        writeln!(file)?;
    }
    console("Data successfully written to Task 1/RK_test.txt");

    // Eigenenergies using shooting method \\
    let (energies0,functions0,calculated_e,sign_f) = shooting_method(&r,n2,0,m,settings.bound_ne);
    let (energies10,functions10,_calculated_e,_sign_f) = shooting_method(&r,n2,10,m,settings.bound_ne);

    let mut file = File::create("Task 1/energies.txt")?;
    writeln!(file, "{:<3}{:>1}{:>20}{:>20}", "#", "n", "E_n (l=0) [a.u.]", "E_n (l=10) [a.u.]")?;
    for i in 0..(energies0.len()).max(energies10.len()) {
        write!(file, "{:4.0}", i)?;
        if i < energies0.len() {write!(file, "{:20.12E}", energies0[i])?}
        else {write!(file, "{:>20}", "")?};
        if i < energies10.len() {write!(file, "{:20.12E}",  energies10[i])?};
        writeln!(file)?;
    }
    console("Data successfully written to Task 1/energies.txt");

    // True -> False to turn off the output of eigenfunctions
    if true {
        let mut file = File::create("Task 1/eigenfunctions0.txt")?;
        write!(file, "{:<3}{:>17}", "#", "R [au]")?;
        for i in 0..energies0.len() {
            write!(file, "{:>20}", format!("psi{} [a.u.]",i))?;
        }
        writeln!(file);
        for i in 0..functions0.ncols() {
            write!(file, "{:20.12E}", r.points[i])?;
            for j in 0..energies0.len() {
                write!(file, "{:20.12E}", functions0[(j, i)])?;
            }
            writeln!(file)?;
        }
        console("Data successfully written to Task 1/eigenfunctions0.txt");

        let mut file = File::create("Task 1/eigenfunctions10.txt")?;
        write!(file, "{:<3}{:>17}", "#", "R [au]")?;
        for i in 0..energies10.len() {
            write!(file, "{:>20}", format!("psi{} [a.u.]",i))?;
        }
        writeln!(file);
        for i in 0..functions10.ncols()  {
            write!(file, "{:20.12E}", r.points[i])?;
            for j in 0..energies10.len() {
                write!(file, "{:20.12E}", functions10[(j, i)])?;
            }
            writeln!(file)?;
        }
        console("Data successfully written to Task 1/eigenfunctions10.txt");
    }

    let mut file = File::create("Test/sign_function.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "e [au]", "f [a.u.]")?;
    for i in 0..calculated_e.len() {
        write!(file, "{:20.12E}", calculated_e[i])?;
        write!(file, "{:20.12E}", sign_f[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to Test/sign_function.txt");

    let energies_control0 = Eigenstates::fourier_dvr(&r, &n2, m,0 ,400, 200)?;
    let energies_control10 = Eigenstates::fourier_dvr(&r, &n2, m,10 ,400, 200)?;
    let mut file = File::create("Task 1/energies_control.txt")?;
    writeln!(file, "{:<3}{:>1}{:>20}{:>20}", "#", "n", "E_n (l=0) [a.u.]", "E_n (l=10) [a.u.]")?;
    for i in 0..200 {
        write!(file, "{:4.0}", i)?;
        if energies_control0.e[i] < 0.0 { write!(file, "{:20.12E}", energies_control0.e[i])? };
        if energies_control10.e[i] < 0.0 { write!(file, "{:20.12E}", energies_control10.e[i])? };
        if energies_control0.e[i] > 0.0 || energies_control10.e[i] > 0.0 {
            break;
        }
        writeln!(file)?;
    }
    console("Data successfully written to Task 1/energies_control.txt");


    // Interpolation output \\
    let n_start = 10;
    let n_final = 2000;
    let n_splines = 8;
    let n_points: Vec<usize>= (0..n_splines).map(|i| {
        let start = (n_start as f64).ln();
        let end = (n_final as f64).ln();
        let t = i as f64 / (n_splines - 1) as f64;
        (start + (end - start) * t).exp().round() as usize
    }).collect();
    let mut splines: Vec<Vec<f64>> = Vec::new();
    let mut energies: Vec<Vec<f64>> = Vec::new();
    for &n in n_points.iter() {
        let x_val = Grid::new_n(settings.rmin,settings.rmax,n).points;
        let y_val = x_val.iter().map(|&x| n2(x)).collect();
        let n2_interpolated = create_natural_cubic_spline(x_val, y_val);
        let (energy,_f,_e,_s) = shooting_method(&r,n2_interpolated.clone(),0,m,settings.bound_ne);
        let spline = r.points.iter().map(|&r| n2_interpolated(r)).collect();
        splines.push(spline);
        let energy_vec: Vec<f64> = energy.data.as_slice().to_vec();
        energies.push(energy_vec);
    }

    let mut file = File::create("Task 3/interpolation.txt")?;
    let mut file_e = File::create("Task 3/interpolation_energies.txt")?;
    let mut file_dif = File::create("Task 3/interpolation_energies_difference.txt")?;
    let mut file_dif0 = File::create("Task 3/interpolation_ground_state_difference.txt")?;
    write!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [eV]")?;
    write!(file_e, "{:<3}{:>1}{:>20}", "#", "n", "E_n")?;
    write!(file_dif, "{:<3}{:>1}", "#", "n")?;
    writeln!(file_dif0, "{:<3}{:>17}{:>20}", "#", "n","ΔE_0")?;
    for n in n_points.iter() {
        write!(file, "{:>20}", &format!("s{}(R) [a.u.]",n))?;
        write!(file_e, "{:>20}", &format!("E_n s{} [a.u.]",n))?;
        write!(file_dif, "{:>20}", &format!("ΔE_n s{} [a.u.]",n))?;
    }
    writeln!(file)?;
    writeln!(file_e)?;
    writeln!(file_dif)?;

    for (j,&x) in r.points.iter().enumerate() {
        write!(file, "{:20.12E}", x) ?;
        write!(file, "{:20.12E}", n2(x) * PHYS_H0) ?;
        for i in 0..n_splines {
            write!(file, "{:20.12E}", splines[i][j] * PHYS_H0) ?;
        }
        writeln!(file) ?;
    }
    console("Data successfully written to Task 3/interpolation.txt");

    for j in 0..energies0.len() {
        write!(file_e, "{:4.0}", j)?;
        write!(file_dif, "{:4.0}", j)?;
        write!(file_e, "{:20.12E}", energies0[j])?;
        for i in 0..n_splines {
            if j < energies[i].len() {
                write!(file_e, "{:20.12E}", energies[i][j])?;
                write!(file_dif, "{:20.12E}", (energies[i][j]-energies0[j]).abs())?;
            }
            else {write!(file_e, "{:>20}", "NaN")?;write!(file_dif, "{:>20}", "NaN")?};
        }
        writeln!(file_e)?;
        writeln!(file_dif)?;
    }
    console("Data successfully written to Task 3/interpolation_energies.txt");
    console("Data successfully written to Task 3/interpolation_energies_difference.txt");

    for i in 0..n_splines {
        write!(file_dif0, "{:20.0}", n_points[i])?;
        write!(file_dif0, "{:20.12E}", (energies[i][0]-energies0[0]).abs())?;
        writeln!(file_dif0)?;
    }
    console("Data successfully written to Task 3/interpolation_ground_state_difference.txt");



    // This code was used for testing purposes and is currently commented out.

    /*
    let mut file = File::create("int_test.txt")?;
    for i in 1..19 {
        let mut n = 10;
        for _ in 0..i {n *= 2}
        let step = PI / (n-1) as f64;
        let values: Vec<f64> = (0..n)
            .map(|x| x as f64 * step)
            .map(|x| x.sin())
            .collect();
        let vec: DVector<f64>= DVector::from_vec(values);
        write!(file, "{:20.0}", n)?;
        write!(file, "{:20.12E}", (vec.integrate(step,Trapezoid) -2.0).abs())?;
        write!(file, "{:20.12E}", (vec.integrate(step,Simpson)-2.0).abs())?;
        write!(file, "{:20.12E}", (vec.integrate(step,INT4)-2.0).abs())?;
        writeln!(file)?;
    }

    let mut file = File::create("prim_test.txt")?;
    for i in 1..15 {
        let mut n = 10;
        for _ in 0..i { n *= 2 }
        let step = PI / (n - 1) as f64;

        let values: Vec<f64> = (0..n)
            .map(|x| x as f64 * step)
            .map(|x| x.cos())
            .collect();
        let vec: DVector<f64> = DVector::from_vec(values);

        let sin_values: Vec<f64> = (0..n)
            .map(|x| x as f64 * step)
            .map(|x| x.sin())
            .collect();
        let sin_vec: DVector<f64> = DVector::from_vec(sin_values);

        let diff_vec2 = vec.primitive(step,P2) - &sin_vec;
        let diff_vec4 = vec.primitive(step,P4) - &sin_vec;
        let mean_abs_diff2: f64 = diff_vec2.iter().map(|x| x.abs()).sum::<f64>() / diff_vec2.len() as f64;
        let mean_abs_diff4: f64 = diff_vec4.iter().map(|x| x.abs()).sum::<f64>() / diff_vec4.len() as f64;
        write!(file, "{:20.0}", n)?;
        write!(file, "{:20.12E}", mean_abs_diff2)?;
        write!(file, "{:20.12E}", mean_abs_diff4)?;
        writeln!(file)?;
    }

    let mut file = File::create("complex_prim.txt")?;
    let x = Grid::new_d(0.0,PI,0.01);
    let vec: Vec<_> = x.points.iter().map(|&x|Complex::new(0.0,x).exp()).collect();
    let prim = DVector::from_vec(vec.clone()).primitive(x.d,P4);
    for (i,x) in x.points.iter().enumerate() {
        write!(file, "{:20.12E}", x)?;
        write!(file, "{:20.12E}", vec[i].real())?;
        write!(file, "{:20.12E}", vec[i].imaginary())?;
        write!(file, "{:20.12E}", prim[i].real())?;
        write!(file, "{:20.12E}", prim[i].imaginary())?;
        writeln!(file)?;
    }
     */

    /*
    let x = Grid::new_d(0.0,20.0,0.01);
    let lho_ode = ODE::new(2, |state: &[f64],x| {-state[0]});
    let ini = vec![1.0,0.0];
    let solution = runge_kutta(
        &*lho_ode,
        &ini,
        &x,
        RKMethod::RK4
    );
    let mut file = File::create("RK4.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "u [au]")?;
    for i in 0..x.n {
        write!(file, "{:20.12E}", x.points[i])?;
        write!(file, "{:20.12E}", solution[0][i])?;
        write!(file, "{:20.12E}", solution[1][i])?;
        writeln!(file)?;
    }
    console("Data successfully written to solution.txt");

    let lho_complex = ODE::new(2, |state: &[Complex<f64>], x| {-state[0]});
    let ini_complex = vec![Complex::new(1.0,0.0),Complex::new(0.0,1.0)];
    let solution_complex = runge_kutta(
        &*lho_complex,
        &ini_complex,
        &x,
        RKMethod::RK4
    );
    let mut file = File::create("RK4_complex.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "u [au]")?;
    for i in 0..x.n {
        write!(file, "{:20.12E}", x.points[i])?;
        write!(file, "{:20.12E}", solution_complex[0][i].real())?;
        write!(file, "{:20.12E}", solution_complex[0][i].imaginary())?;
        writeln!(file)?;
    }
    console("Data successfully written to RK4_complex.txt");

    let x = Grid::new_n(0.0,20.0,2000);
    let m = 1.0;

    let v = |r:f64| -> f64 {
        let a = 1.0;
        let v0 = 3.0;
        - v0 * (-r/a).exp()
    };


    let (energies,functions,calculated_e,sign_f) = shooting_method(&x,v,0,m,1000);

    let mut file = File::create("Test/energies.txt")?;
    writeln!(file, "{:<3}{:>1}{:>20}", "#", "n", "E_n [a.u.]")?;
    for i in 0..energies.len() {
        write!(file, "{:4.0}", i)?;
        write!(file, "{:20.12E}", energies[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to energies.txt");


    let mut file = File::create("Test/eigenfunctions.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}{:>20}{:>20}{:>20}{:>20}", "#", "R [au]", "psi_0 [a.u.]", "psi_1 [a.u.]", "psi_2 [a.u.]", "psi_3 [a.u.]", "psi_... [a.u.]")?;
    for i in 0..functions.ncols() {
        write!(file, "{:20.12E}", x.points[i])?;
        for j in 0..energies.len() {
            write!(file, "{:20.12E}", functions[(j,i)])?;
        }
        writeln!(file)?;
    }
    console("Data successfully written to eigenfunctions.txt");

    let mut file = File::create("Test/sign_function.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "e [au]", "f [a.u.]")?;
    for i in 0..calculated_e.len() {
        write!(file, "{:20.12E}", calculated_e[i])?;
        write!(file, "{:20.12E}", sign_f[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to sign_function.txt");

    let psi = DVector::from_vec(calculate_wavefunction(&x, v.clone(), 0, m, -0.411));
    let mut file = File::create("Test/psi_test.txt")?;
    for i in 0..x.n {
        write!(file, "{:20.12E}", x.points[i])?;
        write!(file, "{:20.12E}", psi[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to psi_test.txt");

     */


    Ok(())
}
