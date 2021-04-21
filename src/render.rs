

use glium::{glutin, Surface, implement_vertex, uniform};
use glium::index::PrimitiveType;


#[derive(Copy, Clone)]
pub struct Vertex {
    pub position: [f32; 2],
    pub color: [f32; 3],
}
implement_vertex!(Vertex, position, color);


#[derive(Copy, Clone)]
pub struct InstanceVertex {
    pub world_position: (f32, f32, f32),
    pub size: f32
}

implement_vertex!(InstanceVertex, world_position, size);
pub struct GLObject<T: Copy>
{
    pub vbo: glium::VertexBuffer<T>,
    pub ibo: Option<glium::IndexBuffer<u32>>
}

impl GLObject<Vertex>
{
    pub fn new(display: &glium::Display, vertices: &[Vertex], indices: Option<Vec<u32>>) -> Self
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


pub fn circle(display: &glium::Display, radius: f32, color: [f32; 3], n_points: u32) -> GLObject<Vertex>
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


pub struct Renderer
{
    pub display: glium::Display,
    frame: Option<glium::Frame>
}

#[allow(dead_code)]
impl Renderer
{
    pub fn new(event_loop: &glutin::event_loop::EventLoop<()>) -> Self
    {
        let wb = glutin::window::WindowBuilder::new();
        let cb = glutin::ContextBuilder::new();
        let display = glium::Display::new(wb, cb, &event_loop).unwrap();
        let frame = None;
        
        Self{display, frame}
    }

    pub fn init(&mut self)
    {
        self.frame = Some(self.display.draw());
        self.frame.as_mut().unwrap().clear_color(0.0, 0.0, 0.0, 0.0);
    }

    pub fn swap(&mut self)
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

    pub fn draw_scene<T: Copy>(&mut self,
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


pub struct Scene(pub Vec<Primitive>);

#[allow(dead_code)]
pub enum Primitive
{
    Disk{x: f32, y: f32, r: f32},
    Line{x1: f32, y1: f32, x2: f32, y2: f32}
}

