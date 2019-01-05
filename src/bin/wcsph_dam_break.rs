// external crate imports
use indicatif::{ProgressBar, ProgressStyle};

// crates imports
use prestige::{
    contact_search::{stash_2d, WorldBounds, NNPS},
    physics::sph::{
        kernel::{CubicKernel},
        wcsph::{
            equations::{equation_of_state, continuity_and_momentum_equation, make_accelerations_zero,
                        apply_gravity},
            WCSPH,
        },
    },
    rk2_initialize, rk2_stage_1, rk2_stage_2, WriteOutput,
};
use simple_shapes::{grid_arange, tank};

// std imports
use std::fs;

fn create_entites(spacing: f32) -> (WCSPH, WCSPH) {
    let (xt, yt) = tank(-1., 3., spacing, -1., 4., spacing, 2);
    let (xf, yf) = grid_arange(
        1.03,
        2.17 + spacing / 2.,
        spacing,
        0.0,
        1.0 + spacing / 2.,
        spacing,
    );

    let fluid_particle_no = xf.len();
    let tank_particle_no = xt.len();

    // define density of the particle of the fluid
    let rho_f = 1000.;
    let fluid_nnps_idx = 0;
    let m_single_par_f = rho_f * spacing.powf(2.);

    // create fluid particle array and setup its properties
    let mut fluid = WCSPH::new_with_xyz(xf.clone(), yf.clone(), vec![0.; fluid_particle_no],
                                        fluid_nnps_idx);
    fluid.m = vec![m_single_par_f; fluid_particle_no];
    fluid.h = vec![1.2 * spacing; fluid_particle_no];

    // define density of the particle of the tank
    let rho_t = 1000.;
    let tank_nnps_idx = 1;
    let m_single_par_t = rho_t * spacing.powf(2.);

    // create tank particle array and setup its properties
    let mut tank = WCSPH::new_with_xyz(xt.clone(), yt.clone(), vec![0.; tank_particle_no],
                                       tank_nnps_idx);
    tank.m = vec![m_single_par_t; tank_particle_no];

    (fluid, tank)
}


fn main() {
    let dim = 2;
    // create particles
    let spacing = 0.03;
    let (mut fluid, tank) = create_entites(spacing);

    println!(
        "Fluid particles: {}, tank particles: {}, Total particles: {}",
        fluid.x.len(),
        tank.x.len(),
        fluid.x.len() + tank.x.len()
    );

    // setup nnps
    let world_bounds = WorldBounds::new(-1.1, 3.1, -1.1, 4.1, 0.0, 0.0, spacing);
    let mut nnps = NNPS::new(2, &world_bounds, dim);

    // solver data
    let dt = 1e-4;
    let mut t = 0.;
    let tf = 1.;
    let mut step_no = 0;
    let pfreq = 100;

    let version = env!("CARGO_MANIFEST_DIR");
    let dir_name = version.to_owned() + "/wcsph_dam_break_output";
    let _tmp = fs::create_dir(&dir_name);

    // create a progressbar
    let total_steps = (tf / dt) as u64;
    let pb = ProgressBar::new(total_steps);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] ({eta})")
            .progress_chars("#>-"));

    // get the kernel
    let cubickernel = CubicKernel::new(2).unwrap();
    let speed_of_sound = 10.59;
    let alpha = 0.27;

    while t < tf {
        // stash the particles into the world's cells
        stash_2d(vec![&fluid, &tank], &mut nnps);

        rk2_initialize(&mut vec![&mut fluid]);

        // stage 1
        make_accelerations_zero(&mut fluid.au, &mut fluid.av, &mut fluid.aw, &mut fluid.arho);
        apply_gravity(&mut fluid.au, &mut fluid.av, &mut fluid.aw, 0.0, -9.81, 0.);

        equation_of_state(&mut fluid.p, &fluid.rho, 1000., 7., speed_of_sound);

        continuity_and_momentum_equation(&fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v, &fluid.w,
                                         &fluid.h, &fluid.p, &fluid.rho, &fluid.c,
                                         &mut fluid.arho, &mut fluid.au, &mut fluid.av, &mut fluid.aw,

                                         &tank.x, &tank.y, &tank.z, &tank.u, &tank.v, &tank.w,
                                         &tank.h, &tank.m, &tank.p, &tank.rho, &tank.c, tank.nnps_idx,
                                         alpha, &nnps, &cubickernel);

        continuity_and_momentum_equation(&fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v, &fluid.w,
                                         &fluid.h, &fluid.p, &fluid.rho, &fluid.c,
                                         &mut fluid.arho, &mut fluid.au, &mut fluid.av, &mut fluid.aw,

                                         &fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v,
                                         &fluid.w, &fluid.h, &fluid.m, &fluid.p, &fluid.rho,
                                         &fluid.c, fluid.nnps_idx, alpha, &nnps, &cubickernel);
        rk2_stage_1(&mut vec![&mut fluid], dt);


        // stage 2
        make_accelerations_zero(&mut fluid.au, &mut fluid.av, &mut fluid.aw, &mut fluid.arho);
        apply_gravity(&mut fluid.au, &mut fluid.av, &mut fluid.aw, 0.0, -9.81, 0.);

        equation_of_state(&mut fluid.p, &fluid.rho, 1000., 7., speed_of_sound);

        continuity_and_momentum_equation(&fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v, &fluid.w,
                                         &fluid.h, &fluid.p, &fluid.rho, &fluid.c,
                                         &mut fluid.arho, &mut fluid.au, &mut fluid.av, &mut fluid.aw,

                                         &tank.x, &tank.y, &tank.z, &tank.u, &tank.v, &tank.w,
                                         &tank.h, &tank.m, &tank.p, &tank.rho, &tank.c, tank.nnps_idx,
                                         alpha, &nnps, &cubickernel);

        continuity_and_momentum_equation(&fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v, &fluid.w,
                                         &fluid.h, &fluid.p, &fluid.rho, &fluid.c,
                                         &mut fluid.arho, &mut fluid.au, &mut fluid.av, &mut fluid.aw,

                                         &fluid.x, &fluid.y, &fluid.z, &fluid.u, &fluid.v,
                                         &fluid.w, &fluid.h, &fluid.m, &fluid.p, &fluid.rho,
                                         &fluid.c, fluid.nnps_idx, alpha, &nnps, &cubickernel);
        rk2_stage_2(&mut vec![&mut fluid], dt);

        if step_no % pfreq == 0 {
            // println!("{}", step_no);
            let filename = format!("{}/tank_{}.vtk", &dir_name, step_no);
            tank.write_vtk(filename);

            let filename = format!("{}/fluid_{}.vtk", &dir_name, step_no);
            fluid.write_vtk(filename);
            // println!("{} ", step_no);
        }
        step_no += 1;
        t += dt;

        // progressbar increment
        pb.inc(1);
    }
    pb.finish_with_message("Simulation succesfully completed");
}
