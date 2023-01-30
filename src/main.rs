//use std::slice::range;

use fltk::image::PngImage;
use fltk::{app, 
    button::Button, 
    frame::Frame, 
    prelude::*, 
    window::Window, 
    input::FloatInput, 
    group::Flex};
use plotters::prelude::*;
use plotters_bitmap::BitMapBackend;
use std::error::Error;

const W: usize = 400;
const H: usize = 200;
//const cells: usize = x_size*y_size*3;
use std::fs;


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

fn decimate_vector(xs: &Vec<f64>, target_len: usize) -> Vec<f64> {
    let current_count = xs.len();
    let ratio = (current_count+1) / target_len;
    let indices_to_pull: Vec<usize> = (0..target_len).into_iter().map(|x| x * ratio).collect();
    indices_to_pull.into_iter().map(|i| xs[i]).collect()
}

#[test]
fn test_decimate_vector() {
    let test_vec = vec![1.0,2.0,3.0,4.0];
    let small_vec = decimate_vector(&test_vec, 2);
    assert_eq!(small_vec.len(), 2);
    assert_eq!(small_vec, vec![1.0,3.0])
}

fn decimate_temps(temps: &(Vec<f64>, Vec<f64>)) -> (Vec<f64>, Vec<f64>){
    (decimate_vector(&temps.0, W/2), 
    decimate_vector(&temps.1, W/2))
}

fn run(t_initial: &f64, thickness: &f64, subdiv_count:&f64, k: &f64, specific_heat: &f64, density: &f64, t_fluid: &f64, h: &f64, duration: &f64, dt: &f64) -> (Vec<f64>,Vec<f64>){
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
    (inner_temps, outer_temps)
}

fn make_chart(inner_temps: Option<Vec<f64>>, outer_temps:Option<Vec<f64>>, timestep:&f64){

    let drawing_area = BitMapBackend::new("tmp.png", (W as u32, H as u32)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();
    let mut chart_builder = ChartBuilder::on(&drawing_area);
    chart_builder.margin(10)
    .x_label_area_size(30)
    .y_label_area_size(50);
    let mut chart_context = match &inner_temps {
        Some(temps) => {
            chart_builder.build_cartesian_2d(
                0.0..temps.len() as f64 * timestep, 
                0.0..1.05*temps[0]).unwrap()       
        }, 
        None => {
            chart_builder.build_cartesian_2d(
                0.0..60.0, 
                0.0..200.0).unwrap()        
        }
    };
    chart_context.configure_mesh()
    .y_desc("Temperature [C]")
    .x_desc("Time [seconds]").draw().unwrap();
    
    
    if let Some(temps) = inner_temps{ 
        chart_context.draw_series(LineSeries::new(
            (0..).zip(
                temps.iter()).map(
                    |(time, temp)|   (time as f64 * timestep, temp.clone())), RED)).unwrap()
                    .label("Inner Wall")
                    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }

    if let Some(temps) =  outer_temps {
        chart_context.draw_series(LineSeries::new(
            (0..).zip(
                temps.iter()).map(
                    |(time, temp)|   (time as f64 * timestep, temp.clone())), BLACK)).unwrap()
                    .label("Outer Wall")
                    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));    
        chart_context.configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .unwrap();
    }   

}
 
fn float_input_row(label: &String, pad: &i32) -> FloatInput {
    let mut row = Flex::default_fill().row();
    row.set_pad(pad.clone());
    Frame::default().with_label(&label);

    let thickness_input = FloatInput::default();
    row.end();
    thickness_input
}

fn delete_temp_file() -> std::io::Result<()> {
    fs::remove_file("tmp.png")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>>{
    //let mut buf = vec![0u8; W * H * 3];
    //let _ = make_blank_chart(&mut image_data);
    let app = app::App::default();
    let mut wind = Window::default().with_size(500, 800).with_label("Quenching Temperatures");
    wind.make_resizable(true);

    let pad = 5;
    let mut main_column = Flex::default_fill()
        .column();
    main_column.set_margin(pad.clone());
    main_column.set_pad(pad.clone());

    // Geometry label
    Frame::default().with_label("Geometry");
    
    // Thickness row
    let mut thickness_input = float_input_row(&"Material Thickness [mm]".to_string(), &pad);
    thickness_input.set_value(&"28.0");

    // physics label
    Frame::default().with_label("Physics");

    // Initial Temp row
    let mut initial_temp_input = float_input_row(&"Initial Material Temperature [C]".to_string(), &pad); 
    initial_temp_input.set_value(&"280.0");

    // Fluid temp row
    let mut fluid_temp_input = float_input_row(&"Fluid Temperature [C]".to_string(), &pad);
    fluid_temp_input.set_value(&"10.0");

    // Convective coefficient row
    let mut convective_input = float_input_row(&"Convective Coefficient [W/m^2*K]".to_string(), &pad);     
    convective_input.set_value(&"2000.0");

    Frame::default().with_label("Simulation Settings");

    // Duration row
    let mut duration_input = float_input_row(&"Duration [Seconds]".to_string(), &pad);
    duration_input.set_value(&"90.0");

    // timestep row
    let mut timestep_input = float_input_row(&"Time Step [Seconds]".to_string(), &pad);
    timestep_input.set_value(&"0.0010");

    // Material Cells Row
    let mut material_subdivs_input = float_input_row(&"Material Subdivisions".to_string(), &pad);
    material_subdivs_input.set_value(&"20");

    let mut button_start = Button::default()
        .with_label("Run Simulation");

    let result_row = Flex::default_fill().row();

    let mut inner_temp_label = Frame::default().with_label("");
    let mut outer_temp_label = Frame::default().with_label("");
    result_row.end();
    
    //let mut image_data: [u8; W*H*3] = [0; W*H*3];
    let _ = make_chart(None, None, &1.0); 
    //println!("make_blanke_chart result: {:?}", res);
    let chart = PngImage::load("tmp.png").expect("Failed to load PNG from disk");
    let mut chart_frame = Frame::default();
    chart_frame.set_image(Some(chart));
    //println!("{:?}",chart_res);

    main_column.set_size(&chart_frame, H as i32);

    main_column.end();

    wind.set_callback(|this_window| {
        let _=delete_temp_file();
        this_window.hide();
    });
    wind.end();
    wind.show();

    // Simulation Params
    let material_k = 40.; // W/m*K
    let specific_heat = 500.; // J/kg*K
    let density = 7850.; // kg/m^3

    /* Event handling */
    button_start.set_callback(move |_| {
        let temps = run(&parse_float_input(&initial_temp_input),
            &(parse_float_input(&thickness_input)/1000.),
            &parse_float_input(&material_subdivs_input),
            &material_k,
            &specific_heat,
            &density,
            &parse_float_input(&fluid_temp_input),
            &parse_float_input(&convective_input),
            &parse_float_input(&duration_input),
            &parse_float_input(&timestep_input));
        
        let inner_temps = &temps.0;
        let outer_temps = &temps.1;
        let inner_temp_string = format!("Inner Wall Temperature: {:.3}",inner_temps[inner_temps.len()-1]);
        inner_temp_label.set_label(&inner_temp_string);
        let outer_temp_string = format!("Outer Wall Temperature: {:.3}",outer_temps[outer_temps.len()-1]);
        outer_temp_label.set_label(&outer_temp_string);

        let sparse_temps = decimate_temps(&temps);
        let sparse_timestep = &parse_float_input(&duration_input) /(W/2) as f64;
        let _ = make_chart(
            Some(sparse_temps.0), 
            Some(sparse_temps.1),
            &sparse_timestep);

        let chart = PngImage::load("tmp.png").expect("Failed to load PNG from disk");
        chart_frame.set_image(Some(chart));
        chart_frame.redraw();    
    });
    app.run().unwrap();
Ok(())
}
