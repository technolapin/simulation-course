
#[macro_use]
extern crate glium;

use glium::{glutin, Surface};
use glium::index::PrimitiveType;

const SIZE_MASS_RATIO: f32 = 0.1;

#[derive(Copy, Clone)]
struct Vertex {
    position: [f32; 2],
    color: [f32; 3],
}
implement_vertex!(Vertex, position, color);


#[derive(Copy, Clone)]
struct InstanceVertex {
    world_position: (f32, f32, f32),
    size: f32
}

implement_vertex!(InstanceVertex, world_position, size);


struct GLObject<T: Copy>
{
    vbo: glium::VertexBuffer<T>,
    ibo: Option<glium::IndexBuffer<u32>>
}

impl GLObject<Vertex>
{
    fn new(display: &glium::Display, vertices: &[Vertex], indices: Option<Vec<u32>>) -> Self
    {
        
        let vbo = glium::VertexBuffer::new(display,
                                           vertices).unwrap();
        
        // building the index buffer
        let ibo = indices.map(|inds| glium::IndexBuffer::new(display,
                                                             PrimitiveType::TrianglesList,
                                                             inds.as_slice()).unwrap());
        
        Self {vbo, ibo}
    }
}


fn circle(display: &glium::Display, radius: f32, color: [f32; 3], n_points: u32) -> GLObject<Vertex>
{
    let pas = std::f32::consts::PI*2.0/(n_points as f32);
    let mut vertices = vec![];
    let mut indices = vec![];
    
    vertices.push(Vertex{position: [0., 0.], color});
    for i in 0..n_points
    {
        let angle = pas*(i as f32);
        let x = angle.cos()*radius;
        let y = angle.sin()*radius;
        vertices.push(Vertex{position: [x, y], color});

    }

    for i in 0..n_points
    {
        indices.push(0);
        indices.push(i+1);
        indices.push((i+1)%n_points +1);

    }
    GLObject::new(display, vertices.as_slice(), Some(indices))
}


struct Renderer
{
    display: glium::Display,
    frame: Option<glium::Frame>
}

impl Renderer
{
    fn new(event_loop: &glutin::event_loop::EventLoop<()>) -> Self
    {
        let wb = glutin::window::WindowBuilder::new();
        let cb = glutin::ContextBuilder::new();
        let display = glium::Display::new(wb, cb, &event_loop).unwrap();
        let frame = None;
        
        Self{display, frame}
    }

    fn init(&mut self)
    {
        self.frame = Some(self.display.draw());
        self.frame.as_mut().unwrap().clear_color(0.0, 0.0, 0.0, 0.0);
    }

    fn swap(&mut self)
    {
        let frame = self.frame.take();
        if let Some(f) = frame
        {
            f.finish().unwrap();
        }
        
    }

    fn draw<T: Copy>(&mut self,
                     obj: &GLObject<T>,
                     program: &glium::Program)
    {
        // building the uniforms
        let x = 0.0;
        let y = 0.1;
        let uniforms = uniform! {
            matrix: [
                [1.0, 0.0, 0.0, x],
                [0.0, 1.0, 0.0, y],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0f32]
            ]
        };

        if self.frame.is_none()
        {
            self.init();
        }

        if let Some(ibo) = &obj.ibo
        {
            self.frame.as_mut().unwrap().draw(&obj.vbo,
                                     ibo,
                                     &program,
                                     &uniforms,
                                     &Default::default()).unwrap();
        }
        else
        {
            self.frame.as_mut().unwrap().draw(&obj.vbo,
                                     &glium::index::NoIndices(glium::index::PrimitiveType::TrianglesList),
                                     &program,
                                     &uniforms,
                                     &Default::default()).unwrap();

        }
    }

    fn draw_instance<T: Copy>(&mut self,
                              obj: &GLObject<T>,
                              instances: &Vec<InstanceVertex>,
                              program: &glium::Program)
    {
        let inst_buffer = glium::vertex::VertexBuffer::dynamic(&self.display, &instances).unwrap();
        // building the uniforms
        let x = 0.0;
        let y = 0.1;
        let uniforms = uniform! {
            matrix: [
                [1.0, 0.0, 0.0, x],
                [0.0, 1.0, 0.0, y],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0f32]
            ]
        };

        if self.frame.is_none()
        {
            self.init();
        }

        let buffer = (&obj.vbo, inst_buffer.per_instance().unwrap());
        
        if let Some(ibo) = &obj.ibo
        {
            self.frame.as_mut().unwrap().draw(buffer,
                                              ibo,
                                              &program,
                                              &uniforms,
                                              &Default::default()).unwrap();
        }
        else
        {
            self.frame.as_mut().unwrap().draw(buffer,
                                              &glium::index::NoIndices(glium::index::PrimitiveType::TrianglesList),
                                              &program,
                                              &uniforms,
                                              &Default::default()).unwrap();

        }
    }

    fn draw_scene<T: Copy>(&mut self,
                           scene: Scene,
                           palette: &Vec<GLObject<T>>,
                           program: &glium::Program)

    {
        let mut instanciation = palette.iter().map(|_| vec![]).collect::<Vec<_>>();
        
        for prim in scene.0.iter()
        {
            use Primitive::*;
            match *prim
                
            {
                Disk{x, y, r} =>
                {
                    instanciation[0].push(
                        InstanceVertex
                        {
                            world_position: (x, y, 0.),
                            size: r
                        }
                    );
                },
                _ => unimplemented!()
            }
            
        }

        for (mesh, instances) in palette.iter().zip(instanciation.iter())
        {
            self.draw_instance(mesh,
                               instances,
                               program)
            
        }
    }
   
}



mod simu
{
    use nalgebra::{Vector3};
    use rayon::prelude::*;

    pub type Float = f64;

    /// Material Point
    #[derive(Default, Debug)]
    pub struct PMat
    {
        pub mass: Float,
        pub pos: Vector3<Float>,
        pub vel: Vector3<Float>,
        pub forces: Vector3<Float>,
        pub fixed: bool
    }
    impl PMat
    {
        pub fn new(mass: Float, x: Float, y: Float) -> Self
        {
            Self
            {
                mass,
                pos: Vector3::new(x, y, 0.),
                .. Default::default()
            }
        }
    }

    #[derive(Debug)]
    pub enum Link
    {
        Spring
        {
            k: Float,
            l: Float
        },
        SpringTreshIn
        {
            k: Float,
            l: Float,
            s: Float
        },
        SpringTreshOut
        {
            k: Float,
            l: Float,
            s: Float
        },
        Brake
        {
            z: Float
        },
        SpringBrake
        {
            k: Float,
            l: Float,
            z: Float
        },
        FrcConst
        {
            frc: Vector3<Float>
        },
        SpringBrakeBreakable
        {
            k: Float,
            l: Float,
            z: Float,
            s: Float,
            active: bool
        }
            
    }

    impl Link
    {
        fn setup(&mut self, a: &PMat, b: &PMat) -> Vector3<Float>
        {
            use Link::*;
            
            match *self
            {
                Spring{k, l} =>
                {
                    let d = b.pos - a.pos;
                    let dir = d.normalize();
                    k*(d-l*dir)
                },
                SpringTreshIn{k, l, s} =>
                {
                    let d = b.pos - a.pos;
                    if d.norm() < s
                    {
                        let dir = d.normalize();
                        k*(d-l*dir)
                    }
                    else
                    {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                },
                SpringTreshOut{k, l, s} =>
                {
                    let d = b.pos - a.pos;
                    if d.norm() > s
                    {
                        let dir = d.normalize();
                        k*(d-l*dir)
                    }
                    else
                    {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                },
                Brake{z} =>
                {
                    z*(b.vel - a.vel)
                },
                SpringBrake{k, l, z} =>
                {
                    let d = b.pos - a.pos;
                    let dir = d.normalize();

                    k*(d-l*dir) + z*(b.vel - a.vel)
                },
                FrcConst{frc} =>
                {
                    frc
                },
                SpringBrakeBreakable{k, l, z, s, mut active} =>
                {
                    if active
                    {
                        let d = b.pos - a.pos;
                        if d.norm() > s
                        {
                            active = false;
                            return Vector3::new(0.0, 0.0, 0.0);
                        }
                        let dir = d.normalize();
                        
                        k*(d-l*dir) + z*(b.vel - a.vel)
                            
                    }
                    else
                    {
                        Vector3::new(0.0, 0.0, 0.0)
                    }
                }
            }
            
        }
    }
    #[derive(Debug)]
    pub struct World
    {
        pub points: Vec<PMat>,
        pub links: Vec<(Link, usize, usize)>
    }

    impl World
    {
        pub fn new() -> Self
        {
            Self
            {
                points: vec![],
                links: vec![],
            }
        }

        pub fn add_point(&mut self, pt: PMat)
        {
            self.points.push(pt);
        }

        pub fn add_link(&mut self, lk: Link, p1: usize, p2: usize)
        {
            self.links.push((lk, p1, p2));
        }

        fn reset_forces(&mut self)
        {
            for pt in self.points.iter_mut()
            {
                pt.forces = Vector3::new(0.0, 0.0, 0.0);
            }
        }
        
        fn setup_links(&mut self)
        {
            self.reset_forces();
            
            for (link, i1, i2) in self.links.iter_mut()
            {
                if let Link::FrcConst{frc} = link
                {
                    let mut p1 = self.points.get_mut(*i1).unwrap();
                    p1.forces += *frc;
                }
                else
                {
                    let frc =
                    {
                        let p1 = self.points.get(*i1).unwrap();
                        let p2 = self.points.get(*i2).unwrap();

                        link.setup(p1, p2)
                    };
                    let mut p1 = self.points.get_mut(*i1).unwrap();
                    p1.forces += frc;
                    let mut p2 = self.points.get_mut(*i2).unwrap();
                    p2.forces -= frc;
                }
                
            }
            
        }

        fn leapfrog(&mut self, h: Float)
        {
            for pt in self.points.iter_mut().filter(|p| !p.fixed)
            {
                pt.vel += h*pt.forces/pt.mass;
                pt.pos += h*pt.vel;
            }
        }
        
        fn eulerexp(&mut self, h: Float)
        {
            for pt in self.points.iter_mut()
            {
                pt.pos += h*pt.vel;
                pt.vel += h*pt.forces/pt.mass;
            }
        } 

        pub fn make_scene(&self) -> crate::Scene
        {

            let primitives = self.points.par_iter().map(|pt|
            {
                let x = pt.pos[0] as f32;
                let y = pt.pos[1] as f32;
                let r = (pt.mass as f32)*crate::SIZE_MASS_RATIO;

                crate::Primitive::Disk{x, y, r}
            }).collect();

            crate::Scene(primitives)
            
        }

        pub fn run(&mut self, step: Float)
        {
            self.setup_links();
            self.leapfrog(step);
        }

        
    }
    
}


use std::collections::HashMap;
pub struct Scene(pub Vec<Primitive>);

pub enum Primitive
{
    Disk{x: f32, y: f32, r: f32},
    Line{x1: f32, y1: f32, x2: f32, y2: f32}
}


use crate::simu::{World, PMat, Link, Float};

// width and height are nb of cells
fn build_2D_grid(world: &mut World, w: usize, h: usize)
{
    let pas = 1.0/ (w.max(h) as Float);
    let cx = w as Float * pas / 2.0;
    let cy = h as Float * pas / 2.0;

    let pos2index = |i: usize, j: usize| i + j*w;
    let index2pos = |index: usize| (index % w, index / w);

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
    use crate::simu::{World, PMat, Link};
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
        build_2D_grid(&mut world, 100, 100);

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

    let mut actual_tps = 0;
    
    // the main loop
    event_loop.run(move |event, _, control_flow| {


        
        if date_prev_tick.elapsed().as_micros() >= 1000000/f_ech
        { 
            world.run(dt);
            date_prev_tick = Instant::now();
            tick_count += 1;



            if tps_probe_date.elapsed().as_micros() >= tps_probe_resolution
            {
                actual_tps = tick_count*1000000 / tps_probe_resolution;
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
