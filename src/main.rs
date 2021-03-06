
#[macro_use]
extern crate glium;

use glium::{glutin, Surface};
use glium::index::PrimitiveType;


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
    use nalgebra::{Vector2};

    type Float = f64;

    /// Material Point
    #[derive(Default, Debug)]
    pub struct PMat
    {
        pub mass: Float,
        pub pos: Vector2<Float>,
        pub vel: Vector2<Float>,
        pub forces: Vector2<Float>,
        pub fixed: bool
    }
    impl PMat
    {
        pub fn new(mass: Float, x: Float, y: Float) -> Self
        {
            Self
            {
                mass,
                pos: Vector2::new(x, y),
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
        Unconditional
        {
            s: Float
        },
        FrcConst
        {
            frc: Vector2<Float>
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
        fn setup(&mut self, a: &PMat, b: &PMat) -> Vector2<Float>
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
                        Vector2::new(0.0, 0.0)
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
                        Vector2::new(0.0, 0.0)
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
                Unconditional{s} =>
                {
                    unimplemented!()
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
                            return Vector2::new(0.0, 0.0);
                        }
                        let dir = d.normalize();
                        
                        k*(d-l*dir) + z*(b.vel - a.vel)
                            
                    }
                    else
                    {
                        Vector2::new(0.0, 0.0)
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
                pt.forces = Vector2::new(0.0, 0.0);
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
            let mut primitives = vec![];

            for pt in self.points.iter()
            {
                let x = pt.pos[0] as f32;
                let y = pt.pos[1] as f32;
                let r = pt.mass as f32;
                primitives.push(crate::Primitive::Disk{x, y, r});
            }

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
struct Scene(Vec<Primitive>);

enum Primitive
{
    Disk{x: f32, y: f32, r: f32},
    Line{x1: f32, y1: f32, x2: f32, y2: f32}
}


fn main() {
    use crate::simu::{World, PMat, Link};
    let mut world = World::new();

    let nb_points = 100;
    let zoom = 0.2;
    let rope_len = 1.0*zoom;
    let mass = 0.1;
    let pas = rope_len/((nb_points-1) as f64);
    world.add_point(PMat::new(mass, -rope_len/2.0, 0.0));

    
    
    for i in 0..nb_points
    {
        world.add_point(PMat::new(mass, ((i+1) as f64)*pas-rope_len/2.0, 0.0));

        
        world.add_link(
            Link::SpringBrake{k: 1.0, l: pas, z: 0.1},
            i,
            i+1
        );

        world.add_link(
            Link::SpringTreshOut{k: 4.0, l: pas*1.5, s: pas*1.5},
            i,
            i+1
        );
        world.add_link(
            Link::SpringTreshIn{k:  10., l: pas*0.5, s: pas*0.5},
            i,
            i+1
        );

    }

    world.points[0].fixed = true;
    world.points[nb_points].fixed = true;
    
    //world.points[5].pos[1] = 0.3;
    
    let event_loop = glutin::event_loop::EventLoop::new();

    let mut renderer = Renderer::new(&event_loop);

    for i in 0..world.points.len()
    {
        world.add_link(
            Link::FrcConst{frc: nalgebra::Vector2::new(0.0, -0.002)},
            i, i
        );
    }

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

    // Here we draw the black background and triangle to the screen using the previously
    // initialised resources.
    //
    // In this case we use a closure for simplicity, however keep in mind that most serious
    // applications should probably use a function that takes the resources as an argument.

    let palette = vec![obj];

    use std::time::Instant;
    
    let mut date_prev_frame = Instant::now();
        

    // the main loop
    event_loop.run(move |event, _, control_flow| {

        if date_prev_frame.elapsed().as_millis() >= 1000/60
        {
            world.run(0.1);
            let scene = world.make_scene();
            renderer.draw_scene(scene, &palette, &program);
            renderer.swap();
            date_prev_frame = Instant::now();

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
                println!("{:?}", event);
                match event
                {
                    glutin::event::DeviceEvent::Key{..} =>
                    {
                        world.points[0].fixed = false;
                    },
                    
                    _ => ()
                };
                glutin::event_loop::ControlFlow::Poll
            }
            _ => glutin::event_loop::ControlFlow::Poll,
        };
    });
}
