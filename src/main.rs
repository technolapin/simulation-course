mod render;
mod simu;

use render::*;
use crate::simu::{World, PMat, Link, Float};

use glium::{glutin, program};






// width and height are nb of cells
fn build_2d_grid(world: &mut World, w: usize, h: usize)
{
    let pas = 1.0/ (w.max(h) as Float);
    let cx = w as Float * pas / 2.0;
    let cy = h as Float * pas / 2.0;

    let pos2index = |i: usize, j: usize| i + j*w;

    let pos2real = |i: usize, j: usize| ((i as Float)*pas - cx, (j as Float)*pas - cy);

    
    let mass = 0.1;


    // points
    for j in 0..h
    {
        for i in 0..w
        {
            let (x, y) = pos2real(i, j);

            world.add_point(
                PMat::new(mass, x, y)
            );
        }
    }


    // structure
    let k = 1000.;
    let l = pas;
    let z = 0.9;

    // structure horizontal
    for j in 0..h
    {
        for i in 0..(w-1)
        {
            let index0 = pos2index(i, j);
            let index1 = pos2index(i+1, j);
            world.add_link(
                Link::SpringBrake{k, l, z},
                index0,
                index1
            );
        }
    }
    // structure vertical
    for j in 0..(h-1)
    {
        for i in 0..w
        {
            let index0 = pos2index(i, j);
            let index1 = pos2index(i, j+1);
            world.add_link(
                Link::SpringBrake{k, l, z},
                index0,
                index1
            );
        }
    }

    // courbure
    let k = 200.;
    let l = pas*2.;
    let z = 0.2;

    // courbure horizontal
    for i in 0..(w-2)
    {
        for j in 0..h
        {
            let index0 = pos2index(i, j);
            let index1 = pos2index(i+2, j);
            world.add_link(
                Link::SpringBrake{k, l, z},
                index0,
                index1 // right
            );
        }
    }
    // courbure vertical
    for i in 0..w
    {
        for j in 0..(h-2)
        {
            let index0 = pos2index(i, j);
            let index1 = pos2index(i, j+2);
            world.add_link(
                Link::SpringBrake{k, l, z},
                index0,
                index1 // right
            );
        }
    }

    // cisaillement
    let k = 200.;
    let l = pas*(2.0f32.sqrt() as Float);
    let z = 0.2;

    for i in 0..(w-1)
    {
        for j in 0..(h-1)
        {
            // a b
            // c d    
            let a = pos2index(i, j);
            let b = pos2index(i+1, j);
            let c = pos2index(i, j+1);
            let d = pos2index(i+1, j+1);
            
            world.add_link(
                Link::SpringBrake{k, l, z},
                a,
                d // upleft to downright
            );
            world.add_link(
                Link::SpringBrake{k, l, z},
                b,
                c // upright to downleft
            );
        }
    }

    for index in 0..(w*h)
    {
        world.add_link(
            Link::FrcConst{frc: nalgebra::Vector3::new(0.0, -1.0*mass, 0.0)},
            index, index
        );
    }
/*
    // points du haut fixes
    for i in 0..w
    {
        let index = pos2index(i, 0);
        world.points[index].fixed = true;

    }
  */  
    // points de gauche fixes
    for j in 0..h
    {
        let index = pos2index(0, j);
        world.points[index].fixed = true;

    }
    
    
}


fn main() {
  
    let mut world = World::new();


    if false
    {
        let nb_points = 1000;
        let zoom = 0.2;
        let rope_len = 1.0*zoom;
        let mass = 0.01;
        let pas = rope_len/((nb_points-1) as f64);
        let y = 0.5;
        world.add_point(PMat::new(mass, -rope_len/2.0, y));

        let g = 1.0;
        let k = 2000f64;
        let l = pas;
        let z = 0.9;
        for i in 0..nb_points
        {
            world.add_point(PMat::new(mass, ((i+1) as f64)*pas-rope_len/2.0, y));

            world.add_link(
                Link::SpringBrake{k, l, z},
                i,
                i+1
            );
        }

        let k = 500f64;
        let l = pas*2.;
        let z = 0.9;

        for i in 0..nb_points-1
        {
            world.add_link(
                Link::SpringBrake{k, l, z},
                i,
                i+2
            );
        }

        
        world.points[0].fixed = true;
        world.points[nb_points].fixed = true;
        
        //world.points[5].pos[1] = 0.3;
        for i in 0..world.points.len()
        {
            world.add_link(
                Link::FrcConst{frc: nalgebra::Vector3::new(0.0, -1.0*g*mass, 0.0)},
                i, i
            );
        }
    }
    else
    {
        build_2d_grid(&mut world, 30, 30);

    }

    
    let event_loop = glutin::event_loop::EventLoop::new();

    let mut renderer = Renderer::new(&event_loop);


    /*
    let vertices = vec![
        Vertex{position: [-0.5, 0.0], color: [1.0, 0.0, 0.0]},
        Vertex{position: [ 0.5, 0.5], color: [0.0, 1.0, 0.0]},
        Vertex{position: [0.5, -0.5], color: [0.0, 0.0, 1.0]},
    ];
    let indices = vec![0u32, 1, 2];
    let obj = GLObject::new(&renderer.display,
                            vertices.as_slice(),
                            Some(indices)
    );
*/
    let obj = circle(&renderer.display, 0.3, [1.0, 0.2, 0.0], 64);
    // compiling shaders and linking them together
    let program = program!(&renderer.display,
        140 => {
            vertex: "
                #version 140
                uniform mat4 matrix;
                in vec2 position;
                in vec3 color;
                in vec3 world_position;
                in float size;

                out vec3 vColor;
                void main() {
                    gl_Position = (vec4(world_position, 0.0) + vec4(position*size, 0.0, 1.0)) * matrix;
                    vColor = color;
                }
            ",

            fragment: "
                #version 140
                in vec3 vColor;
                out vec4 f_color;
                void main() {
                    f_color = vec4(vColor, 1.0);
                }
            "
        },
    ).unwrap();

    let palette = vec![obj];

    use std::time::Instant;
    
    let mut date_prev_tick = Instant::now();


    let f_ech = 10000; //Hz (= tps)
    let fps = 600; //Hz (= fps)
    let dt = 1./(f_ech as f64);

    let mut tick_count = 0;
    let mut tps_probe_date = Instant::now();
    let tps_probe_resolution = 1000000;

    

    
    // the main loop
    event_loop.run(move |event, _, control_flow| {


        
        if date_prev_tick.elapsed().as_micros() >= 1000000/f_ech
        {
            simu::wind(&mut world, nalgebra::Vector3::new(0., 0., 1.), 0.05, 0.2); 
            world.run(dt);
            date_prev_tick = Instant::now();
            tick_count += 1;



            if tps_probe_date.elapsed().as_micros() >= tps_probe_resolution
            {
                let actual_tps = tick_count*1000000 / tps_probe_resolution;
                tps_probe_date = Instant::now();
                tick_count = 0;
                println!("TPS: {}", actual_tps);
            }


            if tick_count % (f_ech/fps) == 0
            {
                let scene = world.make_scene();
                renderer.draw_scene(scene, &palette, &program);
                renderer.swap();
            }
            
        }

        
        *control_flow = match event {
            glutin::event::Event::WindowEvent { event, .. } => match event {
                // Break from the main loop when the window is closed.
                glutin::event::WindowEvent::CloseRequested => glutin::event_loop::ControlFlow::Exit,
                // Redraw the triangle when the window is resized.
                glutin::event::WindowEvent::Resized(..) => {

                    glutin::event_loop::ControlFlow::Poll
                },
                _ => glutin::event_loop::ControlFlow::Poll,
            },
            glutin::event::Event::DeviceEvent{event, ..} =>
            {
                //println!("{:?}", event);
                match event
                {
                    glutin::event::DeviceEvent::Key{..} =>
                    {
                        world.points[0].fixed = !world.points[0].fixed;
                    },
                    
                    _ => ()
                };
                glutin::event_loop::ControlFlow::Poll
            }
            _ => glutin::event_loop::ControlFlow::Poll,
        };
    });
}
