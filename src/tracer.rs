use std::f32;
use canvas::{rgb, Canvas, Point, RGB};

const INF: f32 = f32::MAX;

#[derive(Debug, Clone)]
struct Sphere {
    center: Point,
    radius: f32,
    color: RGB,
    specular: i32,
}

#[derive(Debug, Clone)]
enum Light {
    Ambient { intensity: f32 },
    Point { intensity: f32, position: Point },
    Directional { intensity: f32, direction: Point },
}

#[derive(Debug)]
struct Scene {
    spheres: Vec<Sphere>,
    lights: Vec<Light>,
}

const BACKGROUND_COLOR: RGB = RGB {
    r: 255,
    g: 255,
    b: 255,
};

pub fn draw(canvas: &mut Canvas) {
    let origin = Point::new(0.0, 0.0, 0.0);
    let scene: Scene = Scene {
        spheres: vec![
            // Shiny red sphere
            Sphere {
                center: Point::new(0.0, -1.0, 3.0),
                radius: 1.0,
                color: rgb(255, 0, 0),
                specular: 500,
            },
            // Shiny blue sphere
            Sphere {
                center: Point::new(2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 0, 255),
                specular: 500,
            },
            // Matte green sphere
            Sphere {
                center: Point::new(-2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 255, 0),
                specular: 10,
            },
            // Very shiny yellow sphere
            Sphere {
                center: Point::new(0.0, -5001.0, 0.0),
                radius: 5000.0,
                color: rgb(255, 255, 0),
                specular: 1000,
            },
        ],
        lights: vec![
            Light::Ambient { intensity: 0.2 },
            Light::Point {
                intensity: 0.6,
                position: Point::new(2.0, 1.0, 0.0),
            },
            Light::Directional {
                intensity: 0.2,
                direction: Point::new(1.0, 4.0, 4.0),
            },
        ],
    };
    for x in (-canvas.width / 2)..(canvas.width / 2) {
        for y in (-canvas.height / 2)..(canvas.height / 2) {
            let d = canvas.to_viewport(x as f32, y as f32);
            let color = trace_ray(&scene, origin, d, canvas.depth as f32, INF);
            canvas.put_pixel(x, y, color);
        }
    }
}

fn trace_ray(scene: &Scene, origin: Point, d: Point, t_min: f32, t_max: f32) -> RGB {
    let mut closest_t = INF;
    let mut closest_sphere = None;
    for sphere in scene.spheres.iter() {
        let (t1, t2) = intersect_ray_sphere(origin, d, sphere);
        if t1 >= t_min && t1 <= t_max && t1 < closest_t {
            closest_t = t1;
            closest_sphere = Some(sphere.clone());
        }
        if t2 >= t_min && t2 <= t_max && t2 < closest_t {
            closest_t = t2;
            closest_sphere = Some(sphere.clone());
        }
    }
    if let Some(sphere) = closest_sphere {
        let p = origin + d.mul(closest_t);
        let mut n = p - sphere.center;
        n = n.div(n.length());
        sphere.color.with_intensity(compute_lighting(scene, p, n, d.mul(-1.0), sphere.specular))
    } else {
        BACKGROUND_COLOR
    }
}

fn dot(p1: Point, p2: Point) -> f32 {
    p1.x * p2.x + p1.y * p2.y + p1.z * p2.z
}

fn intersect_ray_sphere(origin: Point, d: Point, sphere: &Sphere) -> (f32, f32) {
    let center = sphere.center;
    let r = sphere.radius;
    let oc = origin - center;

    let k1 = dot(d, d);
    let k2 = 2.0 * dot(oc, d);
    let k3 = dot(oc, oc) - r.powi(2);

    let discriminant: f32 = k2.powi(2) - 4.0 * k1 * k3;
    if discriminant < 0.0 {
        (INF, INF)
    } else {
        (
            ((-k2 + discriminant.sqrt()) / (2.0 * k1)),
            ((-k2 - discriminant.sqrt()) / (2.0 * k1)),
        )
    }
}

fn compute_lighting(scene: &Scene, p: Point, n: Point, v: Point, s: i32) -> f32 {
    let k = scene.lights.iter().fold(0.0, |acc, light| {
        if let &Light::Ambient { intensity } = light {
            acc + intensity
        } else {
            let (v_light, intensity) = match light {
                &Light::Point {
                    position,
                    intensity,
                } => (position - p, intensity),
                &Light::Directional {
                    direction,
                    intensity,
                } => (direction, intensity),
                _ => unreachable!(),
            };
            // Diffuse
            let n_dot_l = dot(n, v_light);
            let diffuse_i = if n_dot_l > 0.0 {
                acc + (intensity * n_dot_l / (n.length() * v_light.length()))
            } else {
                acc
            };
            if s != -1 {
                let r = n.mul(2.0 * dot(n, v_light)) - v_light;
                let r_dot_v = dot(r, v);
                if r_dot_v > 0.0 {
                    diffuse_i + (intensity * (r_dot_v / (r.length() * v.length())).powi(s))
                } else {
                    diffuse_i
                }
            } else {
                diffuse_i
            }
        }
    });
    k
}
