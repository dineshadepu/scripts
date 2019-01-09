// [[file:~/phd/code_phd/scripts/README.org::code:tbb_full][code:tbb_full]]
extern crate prestige;
extern crate simple_shapes;

// crates imports
use prestige::{
    contact_search::{stash_3d, WorldBounds, NNPS},
    physics::rigid_body::{
        equations::{apply_gravity, linear_interparticle_force},
        RB3d,
    },
    EulerIntegrator, WriteOutput, setup_progress_bar,
};

// external crate imports
use simple_shapes::{circle_2d, tank_2d};

// std imports
use std::fs;

fn create_particles_tbb(spacing: f32) -> (Vec<f32>, Vec<f32>,
                                          Vec<f32>, Vec<f32>,
                                         Vec<f32>, Vec<f32>,
                                         Vec<f32>, Vec<f32>,) {
    // fix the spacing between the particles
    let layers = 3;
    // create a tank with 26 cm length, 26 cm height
    let (xt, yt) = tank_2d(0.0, 0.1, spacing, 0.0, 0.1, spacing, layers, true);
    // create a cylinder in 2d (that would be a circle),
    let diameter = 0.01; // in meters
    let (xc1, yc1) = circle_2d((0.05, 0.03), diameter / 2., spacing);
    let (xc2, yc2) = circle_2d((0.03, 0.07), diameter / 2., spacing);
    let (xc3, yc3) = circle_2d((0.05, 0.1), diameter / 2., spacing);
    (xt, yt, xc1, yc1, xc2, yc2, xc3, yc3)
}


fn create_entites(spacing: f32) -> (RB3d, RB3d, RB3d, RB3d){
    let (xt, yt, xc1, yc1, xc2, yc2, xc3, yc3) = create_particles_tbb(spacing);
    // create and setup cylinders
    let mut cylinder1 = RB3d::from_xyr(xc1.clone(), yc1, vec![spacing/2.; xc1.len()]);
    let mut cylinder2 = RB3d::from_xyr(xc2.clone(), yc2, vec![spacing/2.; xc2.len()]);
    let mut cylinder3 = RB3d::from_xyr(xc3.clone(), yc3, vec![spacing/2.; xc3.len()]);
    let cylinder_rho = 2700.;
    let cylinder_m = cylinder_rho * spacing.powf(2.);
    // set the mass
    cylinder1.m = vec![cylinder_m; cylinder1.x.len()];
    cylinder1.nnps_idx = 0;
    cylinder1.initialize();
    cylinder2.m = vec![cylinder_m; cylinder2.x.len()];
    cylinder2.nnps_idx = 1;
    cylinder2.initialize();
    cylinder3.m = vec![cylinder_m; cylinder3.x.len()];
    cylinder3.nnps_idx = 2;
    cylinder3.initialize();

    // create and setup tank
    let mut tank = RB3d::from_xyr(xt.clone(), yt, vec![spacing / 2.; xt.len()]);
    let tank_rho = 1051.;
    let tank_m = tank_rho * spacing.powf(2.);

    // set the mass
    tank.m = vec![tank_m; tank.x.len()];
    tank.nnps_idx = 3;
    tank.initialize();

    (cylinder1, cylinder2, cylinder3, tank)
}


fn print_no_part(pars: Vec<&Vec<f32>>) {
    let mut total_pars = 0;
    for x in pars {
        total_pars += x.len();
    }
    println!("Total particles {}", total_pars);
}

fn main() {
    // The diameter of the cylinder is 1 cm, which is 0.01 m. Let's the spacing be
    // 0.05 cm that would be 5 * 1e-5 m.
    let spacing = 1e-3;
    // dimension
    let dim = 2;

    // particles
    let (mut cylinder1, mut cylinder2, mut cylinder3, tank) = create_entites(spacing);

    let kn = 1e5;

    print_no_part(vec![&cylinder1.x, &cylinder2.x, &cylinder3.x, &tank.x]);

    // setup nnps
    let world_bounds = WorldBounds::new(-0.01, 0.3, -0.01, 0.3, 0.0, 0.0, 2. * spacing);
    let mut nnps = NNPS::new(4, &world_bounds, dim);

    // solver data
    let dt = 1e-4;
    let mut t = 0.;
    let tf = 5.0;
    let mut step_no = 0;
    let pfreq = 100;

    let project_root = env!("CARGO_MANIFEST_DIR");
    let dir_name = project_root.to_owned() + "/tbb_1_output";
    let _p = fs::create_dir(&dir_name);

    // create a progress bar
    let total_steps = (tf / dt) as u64;
    let pb = setup_progress_bar(total_steps);
    while t < tf {
        // stash the particles into the world's cells
        stash_3d(vec![&cylinder1, &cylinder2, &cylinder3, &tank], &mut nnps);

        apply_gravity(
            &cylinder1.m, &mut cylinder1.fx, &mut cylinder1.fy, &mut cylinder1.fz,
            0.0, -9.81, 0.0,
        );
        apply_gravity(
            &cylinder2.m, &mut cylinder2.fx, &mut cylinder2.fy, &mut cylinder2.fz,
            0.0, -9.81, 0.0,
        );
        apply_gravity(
            &cylinder3.m, &mut cylinder3.fx, &mut cylinder3.fy, &mut cylinder3.fz,
            0.0, -9.81, 0.0,
        );

        linear_interparticle_force(
            &cylinder1.x, &cylinder1.y, &cylinder1.z,
            &cylinder1.u, &cylinder1.v, &cylinder1.w, &cylinder1.rad,
            &mut cylinder1.fx, &mut cylinder1.fy, &mut cylinder1.fz,

            &tank.x, &tank.y, &tank.z,
            &tank.u, &tank.v, &tank.w, &tank.rad, tank.nnps_idx,

            &nnps, kn, 5.
        );

        // force on cylinder1 due to cylinder 2
        linear_interparticle_force(
            &cylinder1.x, &cylinder1.y, &cylinder1.z,
            &cylinder1.u, &cylinder1.v, &cylinder1.w, &cylinder1.rad,
            &mut cylinder1.fx, &mut cylinder1.fy, &mut cylinder1.fz,

            &cylinder2.x, &cylinder2.y, &cylinder2.z,
            &cylinder2.u, &cylinder2.v, &cylinder2.w, &cylinder2.rad, cylinder2.nnps_idx,

            &nnps, kn, 5.
        );

        // force on cylinder1 due to cylinder 3
        linear_interparticle_force(
            &cylinder1.x, &cylinder1.y, &cylinder1.z,
            &cylinder1.u, &cylinder1.v, &cylinder1.w, &cylinder1.rad,
            &mut cylinder1.fx, &mut cylinder1.fy, &mut cylinder1.fz,

            &cylinder3.x, &cylinder3.y, &cylinder3.z,
            &cylinder3.u, &cylinder3.v, &cylinder3.w, &cylinder3.rad, cylinder3.nnps_idx,

            &nnps, kn, 5.
        );
        // ------------------------------------------

        // ------------------------------------------
        // force on cylinder2 due to tank
        linear_interparticle_force(
            &cylinder2.x, &cylinder2.y, &cylinder2.z,
            &cylinder2.u, &cylinder2.v, &cylinder2.w, &cylinder2.rad,
            &mut cylinder2.fx, &mut cylinder2.fy, &mut cylinder2.fz,

            &tank.x, &tank.y, &tank.z,
            &tank.u, &tank.v, &tank.w, &tank.rad, tank.nnps_idx,

            &nnps, kn, 5.
        );
        // force on cylinder2 due to cylinder 1
        linear_interparticle_force(
            &cylinder2.x, &cylinder2.y, &cylinder2.z,
            &cylinder2.u, &cylinder2.v, &cylinder2.w, &cylinder2.rad,
            &mut cylinder2.fx, &mut cylinder2.fy, &mut cylinder2.fz,

            &cylinder1.x, &cylinder1.y, &cylinder1.z,
            &cylinder1.u, &cylinder1.v, &cylinder1.w, &cylinder1.rad, cylinder1.nnps_idx,

            &nnps, kn, 5.
        );
        // force on cylinder2 due to cylinder 3
        linear_interparticle_force(
            &cylinder2.x, &cylinder2.y, &cylinder2.z,
            &cylinder2.u, &cylinder2.v, &cylinder2.w, &cylinder2.rad,
            &mut cylinder2.fx, &mut cylinder2.fy, &mut cylinder2.fz,

            &cylinder3.x, &cylinder3.y, &cylinder3.z,
            &cylinder3.u, &cylinder3.v, &cylinder3.w, &cylinder3.rad, cylinder3.nnps_idx,

            &nnps, kn, 5.
        );
        // ------------------------------------------

        // ------------------------------------------
        // force on cylinder3 due to tank
        linear_interparticle_force(
            &cylinder3.x, &cylinder3.y, &cylinder3.z,
            &cylinder3.u, &cylinder3.v, &cylinder3.w, &cylinder3.rad,
            &mut cylinder3.fx, &mut cylinder3.fy, &mut cylinder3.fz,

            &tank.x, &tank.y, &tank.z,
            &tank.u, &tank.v, &tank.w, &tank.rad, tank.nnps_idx,

            &nnps, kn, 5.
        );
        // force on cylinder3 due to cylinder 1
        linear_interparticle_force(
            &cylinder3.x, &cylinder3.y, &cylinder3.z,
            &cylinder3.u, &cylinder3.v, &cylinder3.w, &cylinder3.rad,
            &mut cylinder3.fx, &mut cylinder3.fy, &mut cylinder3.fz,

            &cylinder1.x, &cylinder1.y, &cylinder1.z,
            &cylinder1.u, &cylinder1.v, &cylinder1.w, &cylinder1.rad, cylinder1.nnps_idx,

            &nnps, kn, 5.
        );
        // force on cylinder3 due to cylinder 2
        linear_interparticle_force(
            &cylinder3.x, &cylinder3.y, &cylinder3.z,
            &cylinder3.u, &cylinder3.v, &cylinder3.w, &cylinder3.rad,
            &mut cylinder3.fx, &mut cylinder3.fy, &mut cylinder3.fz,

            &cylinder2.x, &cylinder2.y, &cylinder2.z,
            &cylinder2.u, &cylinder2.v, &cylinder2.w, &cylinder2.rad, cylinder2.nnps_idx,

            &nnps, kn, 5.
        );


        cylinder1.euler_stage_1(dt);
        cylinder2.euler_stage_1(dt);
        cylinder3.euler_stage_1(dt);

        if step_no % pfreq == 0 {
            tank.write_vtk(format!("{}/tank_{}.vtk", &dir_name, step_no));
            cylinder1.write_vtk(format!("{}/cylinder1_{}.vtk", &dir_name, step_no));
            cylinder2.write_vtk(format!("{}/cylinder2_{}.vtk", &dir_name, step_no));
            cylinder3.write_vtk(format!("{}/cylinder3_{}.vtk", &dir_name, step_no));
        }
        step_no += 1;
        t += dt;

        // progress bar increment
        pb.inc(1);
    }
    pb.finish_with_message("Simulation succesfully completed");
}
// code:tbb_full ends here
