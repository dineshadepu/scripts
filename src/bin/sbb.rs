// [[file:~/phd/code_phd/scripts/README.org::code:sbb_full][code:sbb_full]]
extern crate indicatif;
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
use indicatif::{ProgressBar, ProgressStyle};
use simple_shapes::{circle_2d, tank_2d};

// std imports
use std::fs;

fn create_particles_sbb(spacing: f32) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    // fix the spacing between the particles
    let layers = 1;
    // create a tank with 26 cm length, 26 cm height
    let (xt, yt) = tank_2d(0.0, 0.26, spacing, 0.0, 0.26, spacing, layers);
    // create a cylinder in 2d (that would be a circle),
    let diameter = 0.01; // in meters
    let (xc, yc) = circle_2d((0.1, 0.1), diameter / 2., spacing);
    (xt, yt, xc, yc)
}

fn create_entites(spacing: f32) -> (RB3d, RB3d){
    let (xt, yt, xc, yc) = create_particles_sbb(spacing);
    // create and setup cylinder
    let mut cylinder = RB3d::from_xyr(xc.clone(), yc, vec![spacing / 2.; xc.len()]);
    let cylinder_rho = 2700.;
    let cylinder_m = cylinder_rho * spacing.powf(2.);
    // set the mass
    cylinder.m = vec![cylinder_m; cylinder.x.len()];
    cylinder.nnps_idx = 0;
    cylinder.initialize();

    // create and setup tank
    let mut tank = RB3d::from_xyr(xt.clone(), yt, vec![spacing / 2.; xt.len()]);
    let tank_rho = 1051.;
    let tank_m = tank_rho * spacing.powf(2.);

    // set the mass
    tank.m = vec![tank_m; tank.x.len()];
    tank.nnps_idx = 1;
    tank.initialize();

    (cylinder, tank)
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
    let spacing = 5. * 1e-5;
    // dimension
    let dim = 2;

    // particles
    let (mut cylinder, tank) = create_entites(spacing);

    let kn = 1e5;

    print_no_part(vec![&cylinder.x, &tank.x]);

    // setup nnps
    let world_bounds = WorldBounds::new(-0.01, 0.3, -0.01, 0.3, 0.0, 0.0, 2. * spacing);
    let mut nnps = NNPS::new(2, &world_bounds, dim);

    // solver data
    let dt = 1e-4;
    let mut t = 0.;
    let tf = 1.;
    let mut step_no = 0;
    let pfreq = 100;

    let project_root = env!("CARGO_MANIFEST_DIR");
    let dir_name = project_root.to_owned() + "/sbb_1_output";
    let _p = fs::create_dir(&dir_name);

    // create a progress bar
    let total_steps = (tf / dt) as u64;
    let pb = setup_progress_bar(total_steps);
    while t < tf {
        // stash the particles into the world's cells
        stash_3d(vec![&cylinder, &tank], &mut nnps);

        apply_gravity(
            &cylinder.m, &mut cylinder.fx, &mut cylinder.fy, &mut cylinder.fz,
            0.0, -9.81, 0.0,
        );
        linear_interparticle_force(
            &cylinder.x, &cylinder.y, &cylinder.z, &cylinder.u,
            &cylinder.v, &cylinder.w, &cylinder.rad, &mut cylinder.fx,
            &mut cylinder.fy, &mut cylinder.fz,

            &tank.x, &tank.y, &tank.z, &tank.u,
            &tank.v, &tank.w, &tank.rad, tank.nnps_idx,

            &nnps,
            kn,
            5.,
        );

        cylinder.euler_stage_1(dt);

        if step_no % pfreq == 0 {
            tank.write_vtk(format!("{}/tank_{}.vtk", &dir_name, step_no));
            cylinder.write_vtk(format!("{}/cylinder_{}.vtk", &dir_name, step_no));
        }
        step_no += 1;
        t += dt;

        // progress bar increment
        pb.inc(1);
    }
    pb.finish_with_message("Simulation succesfully completed");
}
// code:sbb_full ends here
