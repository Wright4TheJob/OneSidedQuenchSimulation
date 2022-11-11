use fltk::{app, button::Button, frame::Frame, prelude::*, window::Window, input::FloatInput, group::Flex};
fn main() {
    let app = app::App::default();
    let mut wind = Window::default().with_size(500, 800).with_label("Quenching Temperatures");
    let mut flex = Flex::default()
        .with_size(240,800)
        .with_pos(250,0)
        .column();
    let mut geometry_label = Frame::default().with_label("Geometry");
    let mut thickness_input = FloatInput::default()
        .with_label("Material Thickness [mm]");
    thickness_input.set_value(&"28.0");
    let mut pysics_label = Frame::default().with_label("Physics");
    let mut initial_temp_input = FloatInput::default()
        .with_label("Initial Material Temperature [C]");
    initial_temp_input.set_value(&"280.0");
    let mut fluid_temp_input = FloatInput::default()
        .with_label("Fluid Temperature [C]");
    fluid_temp_input.set_value(&"10.0");
    let mut convective_input = FloatInput::default()
        .with_label("Convective Coefficient [W/m^2*K]");
    convective_input.set_value(&"2000.0");
    let mut simulation_label = Frame::default().with_label("Simulation Settings");
    let mut duration_input = FloatInput::default()
        .with_label("Duration [Seconds]");
    duration_input.set_value(&"90.0");
    let mut timestep_input = FloatInput::default()
        .with_label("Time Step [Seconds]");
    timestep_input.set_value(&"0.0010");
    let mut material_subdivs_input = FloatInput::default()
        .with_label("Material Subdivisions");
    material_subdivs_input.set_value(&"20");
    let mut button_start = Button::default()
        .with_label("Run Simulation");
    let mut inner_temp_label = Frame::default().with_label("");
    let mut outer_temp_label = Frame::default().with_label("");
    flex.end();
    wind.end();
    wind.show();

    // Simulation Params
    let mut subdiv_count = 10.;
    let mut thickness = 0.028; // meters
    let mut temp_initial = 250.0; // Celcius
    let mut material_k = 40.; // W/m*K
    let mut specific_heat = 500.; // J/kg*K
    let mut density = 7850.; // kg/m^3
    let mut convective_h = 2000.; // W/m^2*K
    let mut fluid_temp = 20.; // C
    let mut total_duration = 60.; //seconds
    let mut timestep = 0.001; // seconds

    /* Event handling */
    button_start.set_callback(move |_|
        run(&parse_float_input(&initial_temp_input),
            &(parse_float_input(&thickness_input)/1000.),
            &parse_float_input(&material_subdivs_input),
            &material_k,
            &specific_heat,
            &density,
            &parse_float_input(&fluid_temp_input),
            &parse_float_input(&convective_input),
            &parse_float_input(&duration_input),
            &parse_float_input(&timestep_input),
            &mut inner_temp_label,
            &mut outer_temp_label));
    app.run().unwrap();
}

fn parse_float_input(field: &FloatInput) -> f64 {
    let val: f64 = if field.value().is_empty() {
        0.0
    } else {
        field.value().parse().unwrap_or(0.0)
    };
    return val;
}

#[test]
fn outside_test() {
    let result = true;
    assert!(result);
}

#[test]
fn cells_initialized_correct_length() {
    let desired = 10;
    let cells = initialize_cells(10, 200.0);
    assert_eq!(desired, cells.len());
    assert_eq!(200.0, cells[0]);
}

fn heat_flow_between_cells(t1: &f64, t2: &f64, dx: &f64, k: &f64) -> f64 {
    let dt = t2 - t1; // degrees C or K
    return -k*dt/dx; // w/m^2
}

#[test]
fn test_reverse_heat_flow() {
    let heat = heat_flow_between_cells(&10.0, &0.0, &1.0, &1.0);
    assert_eq!(heat, 10.0)
}
#[test]
fn test_basic_heat_flow() {
    let heat = heat_flow_between_cells(&1.0, &2.0, &1.0, &1.0);
    assert_eq!(heat, -1.0)
}

#[test]
fn test_insulating_heat_flow() {
    let heat = heat_flow_between_cells(&1.0, &2.0, &1.0, &0.02);
    assert_eq!(heat, -0.02)
}

#[test]
fn test_thick_material_heat_flow() {
    let heat = heat_flow_between_cells(&1.0, &2.0, &5.0, &0.5);
    assert_eq!(heat, -0.1)
}

fn convective_heat_flow(t_fluid: &f64, t_wall: &f64, h: &f64) -> f64{
    return h*(t_fluid - t_wall)
}

#[test]
fn test_convection_simple() {
    assert_eq!(convective_heat_flow(&1.0, &0.0, &1.0), 1.0)
}

#[test]
fn test_convection_negative() {
    assert_eq!(convective_heat_flow(&0.0, &5.0, &1.0), -5.0)
}

#[test]
fn test_convection_negative_high_coeff() {
    assert_eq!(convective_heat_flow(&0.0, &5.0, &10.0), -50.0)
}

fn new_temperature(t_initial: &f64,
                    specific_heat: &f64,
                    mass: &f64,
                    net_heat_flow: &f64,
                    dt: &f64) -> f64 {
    let thermal_mass = mass*specific_heat; // J/K
    let net_energy = net_heat_flow * dt; // J
    let t_final = t_initial + net_energy/thermal_mass;
    return t_final;
}

#[test]
fn test_change_temp_simple() {
    let dt = new_temperature(&0.0, &1., &1., &10., &1.);
    assert_eq!(dt, 10.);
}

fn get_mass(density: f64, dx: f64) -> f64 {
    return density*dx;
}

fn conductive_flows(ts: &Vec<f64>, dx: &f64, k: &f64) -> Vec<f64> {
    let mut flows = Vec::with_capacity(ts.len()-1);
    for i in 0..ts.len()-1 {
        flows.push(heat_flow_between_cells(&ts[i], &ts[i+1], &dx, &k));
    }
    return flows;
}

#[test]
fn test_conductive_flows_simple() {
    let ts = vec![1., 1.];
    let flows = conductive_flows(&ts, &0.1, &1.0);
    assert_eq!(0., flows[0])
}

fn net_heat_flows(ts:&Vec<f64>, conductive_flows: &Vec<f64>, h: &f64, t_fluid: &f64) -> Vec<f64> {
    let mut net_flows = Vec::with_capacity(ts.len());
    for i in 0..ts.len() {
        let left;
        let mut right = 0.;
        if i == 0 {
            left = convective_heat_flow(&t_fluid, &ts[0], &h);
            right = conductive_flows[0];
        } else if i == ts.len()-1 {
            left = conductive_flows[i-1];
        } else {
            left = conductive_flows[i-1];
            right = conductive_flows[i];
        }
        net_flows.push(left - right);
    }
    return net_flows
}

fn update_temperatures(ts: &Vec<f64>, net_heat: &Vec<f64>, specific_heat: &f64, mass: &f64, dt: &f64) -> Vec<f64> {
    let mut new_temps = Vec::with_capacity(ts.len());
    for i in 0..ts.len() {
        new_temps.push(new_temperature(&ts[i], &specific_heat, &mass, &net_heat[i], &dt));
    }
    return new_temps;
}

fn run(t_initial: &f64, thickness: &f64, subdiv_count:&f64, k: &f64, specific_heat: &f64, density: &f64, t_fluid: &f64, h: &f64, duration: &f64, dt: &f64, inner_label: &mut Frame, outer_label: &mut Frame){
    let dx = thickness.clone()/subdiv_count.clone();
    let n = subdiv_count.clone() as i64;
    let cell_mass = get_mass(density.clone(), dx);
    let mut temps = vec![t_initial.clone(); n.clone() as usize];
    let mut outer_temps: Vec<f64> = Vec::with_capacity(temps.len());
    let mut inner_temps: Vec<f64> = Vec::with_capacity(temps.len());
    let mut t = 0.;

    outer_temps.push(temps[0]);
    inner_temps.push(temps[n.clone() as usize-1]);

    while t <= duration.clone() {
        let conductions = conductive_flows(&temps, &dx, &k);
        let nets = net_heat_flows(&temps, &conductions, &h, &t_fluid);
        let new_temps = update_temperatures(&temps, &nets, &specific_heat, &cell_mass, &dt);
        temps = new_temps;
        outer_temps.push(temps[0].clone());
        inner_temps.push(temps[n.clone() as usize-3].clone());
        t = t + dt;
    }
    let mut inner_temp_string = format!("Inner Wall Temperature: {:.3}",inner_temps[inner_temps.len()-1]);
    inner_label.set_label(&inner_temp_string);
    let mut outer_temp_string = format!("Outer Wall Temperature: {:.3}",outer_temps[outer_temps.len()-1]);
    outer_label.set_label(&outer_temp_string);
}
