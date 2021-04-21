const SIZE_MASS_RATIO: f32 = 0.1;


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

#[allow(dead_code)]
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
        
        match self
        {
            Spring{k, l} =>
            {
                let d = b.pos - a.pos;
                let dir = d.normalize();
                *k * (d - *l * dir)
            },
            SpringTreshIn{k, l, s} =>
            {
                let d = b.pos - a.pos;
                if d.norm() < *s
                {
                    let dir = d.normalize();
                    *k * (d - *l * dir)
                }
                else
                {
                    Vector3::new(0.0, 0.0, 0.0)
                }
            },
            SpringTreshOut{k, l, s} =>
            {
                let d = b.pos - a.pos;
                if d.norm() > *s
                {
                    let dir = d.normalize();
                    *k*(d-*l*dir)
                }
                else
                {
                    Vector3::new(0.0, 0.0, 0.0)
                }
            },
            Brake{z} =>
            {
                *z * (b.vel - a.vel)
            },
            SpringBrake{k, l, z} =>
            {
                let d = b.pos - a.pos;
                let dir = d.normalize();

                *k * (d - *l * dir) + *z * (b.vel - a.vel)
            },
            FrcConst{frc} =>
            {
                *frc
            },
            SpringBrakeBreakable{k, l, z, s, active} =>
            {
                if *active
                {
                    let d = b.pos - a.pos;
                    if d.norm() > *s
                    {
                        *active = false;
                        return Vector3::new(0.0, 0.0, 0.0);
                    }
                    let dir = d.normalize();
                    
                    *k * (d - *l * dir) + *z * (b.vel - a.vel)
                        
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
                let p1 = self.points.get_mut(*i1).unwrap();
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
                let p1 = self.points.get_mut(*i1).unwrap();
                p1.forces += frc;
                let p2 = self.points.get_mut(*i2).unwrap();
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
    #[allow(dead_code)]    
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
                                                        let r = (pt.mass as f32)*SIZE_MASS_RATIO;

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

pub fn wind(world: &mut World, direction: Vector3<Float>, probability: f32, intensity: Float)
{
    for _ in 0..((world.points.len() as f32 * probability) as usize)
    {
        let index = rand::random::<usize>() % world.points.len();
        world.points[index].vel += direction * intensity;
    }
}


