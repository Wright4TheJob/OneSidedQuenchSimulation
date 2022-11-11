fn main() {
    // Simulation Params
    let subdiv_count = 10;
    let thickness = 0.028; // meters
    let temp_initial = 250.0; // Celcius
    let material_k = 40.; // W/m*K
    let specific_heat = 500.; // J/kg*K
    let density = 7850.; // kg/m^3
    let convective_h = 2000.; // W/m^2*K
    let fluid_temp = 20.; // C
    let total_duration = 60.; //seconds
    let timestep = 0.1; // seconds
    let dx = thickness/subdiv_count as f64;
    let cell_mass = mass(density, dx);
    let (_final_temps, outer_temps, inner_temps) = iterate(temp_initial, dx, material_k, specific_heat, cell_mass, fluid_temp, convective_h, subdiv_count, total_duration, timestep);
    println!("Inner Wall Temperature: {:.3}",inner_temps[inner_temps.len()-1]);
    println!("Outer Wall Temperature: {:.3}",outer_temps[outer_temps.len()-1]);
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

fn mass(density: f64, dx: f64) -> f64 {
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

fn iterate(t_initial: f64, dx: f64, k: f64, specific_heat: f64, mass: f64, t_fluid: f64, h: f64, n: i32, duration: f64, dt: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut temps = vec![t_initial; n as usize];
    let mut outer_temps: Vec<f64> = Vec::with_capacity(temps.len());
    let mut inner_temps: Vec<f64> = Vec::with_capacity(temps.len());
    let mut t = 0.;

    outer_temps.push(temps[0]);
    inner_temps.push(temps[n as usize-1]);

    while t <= duration {
        let conductions = conductive_flows(&temps, &dx, &k);
        let nets = net_heat_flows(&temps, &conductions, &h, &t_fluid);
        let new_temps = update_temperatures(&temps, &nets, &specific_heat, &mass, &dt);
        temps = new_temps;
        outer_temps.push(temps[0].clone());
        inner_temps.push(temps[n as usize-3].clone());
        t = t + dt;
    }
    return (temps, outer_temps, inner_temps);
}
