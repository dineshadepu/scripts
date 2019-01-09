use simple_shapes::{grid_arange, hollow_box_2d};

pub struct ShearDrivenCavityProperties {
    pub xf: Vec<f32>,
    pub yf: Vec<f32>,
    pub uf: Vec<f32>,
    pub vf: Vec<f32>,
    pub xb: Vec<f32>,
    pub yb: Vec<f32>,
    pub ub: Vec<f32>,
    pub vb: Vec<f32>,
    pub speed_of_sound: f32,
    pub art_visc: f32,
    pub fluid_side_len: f32,
    pub fluid_spacing: f32,
    pub boundary_spacing: f32,
    pub boundary_layers: usize,
    pub velocity_top_layer: f32,
}

impl ShearDrivenCavityProperties {
    pub fn default() -> Self {
        // spacing of particles
        let fluid_side_len = 0.001;
        let fluid_spacing = 2.5 * 1e-5;
        let speed_of_sound = 0.25;
        let fluid_side_len = 0.001;
        let art_visc = 1.95;
        let boundary_spacing = 2.5 * 1e-5;
        let boundary_layers = 2;

        // create a grid with side length of 0.001m
        let (xf, yf) = grid_arange(
            0.,
            fluid_side_len,
            fluid_spacing,
            0.,
            fluid_side_len,
            fluid_spacing,
        );

        // create a box which closes the fluid.
        let (xb, yb) = hollow_box_2d(
            0.,
            0.001,
            fluid_spacing,
            0.,
            0.001,
            fluid_spacing,
            1,
            true,
        );

        ShearDrivenCavityProperties {
            xb: xb.clone(),
            yb: yb.clone(),
            xf: xf.clone(),
            yf: yf.clone(),
            uf: vec![0.; xf.len()],
            vf: vec![0.; xf.len()],
            ub: vec![0.; xb.len()],
            vb: vec![0.; xb.len()],
            speed_of_sound: speed_of_sound,
            art_visc: art_visc,
            fluid_side_len: fluid_side_len,
            fluid_spacing: fluid_spacing,
            boundary_spacing: boundary_spacing,
            boundary_layers: boundary_layers,
            velocity_top_layer: 0.001,
        }
    }
    pub fn get_positions(&self) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>){
        (self.xf.clone(), self.yf.clone(), self.xb.clone(), self.xb.clone())
    }
}
