use std::f32;
use canvas::{rgb, Canvas, Point, RGB};

const INF: f32 = f32::MAX;

#[derive(Debug, Clone)]
struct Sphere {
    center: Point,
    radius: f32,
    color: RGB,
    specular: i32,
    reflective: f32,
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
    r: 0,
    g: 0,
    b: 0,
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
                reflective: 0.2,
            },
            // Shiny blue sphere
            Sphere {
                center: Point::new(2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 0, 255),
                specular: 500,
                reflective: 0.3,
            },
            // Matte green sphere
            Sphere {
                center: Point::new(-2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 255, 0),
                specular: 10,
                reflective: 0.4,
            },
            // Very shiny yellow sphere
            Sphere {
                center: Point::new(0.0, -5001.0, 0.0),
                radius: 5000.0,
                color: rgb(255, 255, 0),
                specular: 1000,
                reflective: 0.5,
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
            let color = trace_ray(&scene, origin, d, canvas.depth as f32, INF, 3);
            canvas.put_pixel(x, y, color);
        }
    }
}

fn closest_intersection(
    scene: &Scene,
    origin: Point,
    d: Point,
    t_min: f32,
    t_max: f32,
) -> (Option<Sphere>, f32) {
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
    (closest_sphere, closest_t)
}

fn trace_ray(scene: &Scene, origin: Point, d: Point, t_min: f32, t_max: f32, depth: u32) -> RGB {
    let (closest_sphere, closest_t) = closest_intersection(scene, origin, d, t_min, t_max);
    if closest_sphere.is_none() {
        return BACKGROUND_COLOR;
    }
    let sphere = closest_sphere.unwrap();

    // Compute local color
    let p = origin + d.mul(closest_t);
    let mut n = p - sphere.center;
    n = n.div(n.length());
    let local_color = sphere.color.with_intensity(compute_lighting(scene, p, n, d.mul(-1.0), sphere.specular));

    // If we hit the recursion limit or the object is not relective, we're done
    let r = sphere.reflective;
    if depth <= 0 || r <= 0.0 {
        return local_color;
    }

    let ray = reflect_ray(d.mul(-1.0), n);
    let reflected_color = trace_ray(scene, p, ray, 0.1, INF, depth - 1);

    local_color.with_intensity(1.0 - r) + reflected_color.with_intensity(r)
}

fn reflect_ray(ray: Point, n: Point) -> Point {
    n.mul(2.0 * n.dot(ray)) - ray
}

fn intersect_ray_sphere(origin: Point, d: Point, sphere: &Sphere) -> (f32, f32) {
    let center = sphere.center;
    let r = sphere.radius;
    let oc = origin - center;

    let k1 = d.dot(d);
    let k2 = 2.0 * oc.dot(d);
    let k3 = oc.dot(oc) - r.powi(2);

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
    scene.lights.iter().fold(0.0, |mut total_intensity, light| {
        if let &Light::Ambient { intensity } = light {
            return total_intensity + intensity;
        }

        let (v_light, intensity, t_max) = match light {
            &Light::Point {
                position,
                intensity,
            } => (position - p, intensity, 1.0),
            &Light::Directional {
                direction,
                intensity,
            } => (direction, intensity, INF),
            _ => unreachable!(),
        };

        // Shadow check
        let (shadow_sphere, _) = closest_intersection(scene, p, v_light, 0.001, t_max);
        if !shadow_sphere.is_none() {
            return total_intensity;
        }

        // Diffuse
        let n_dot_l = n.dot(v_light);
        if n_dot_l > 0.0 {
            total_intensity += intensity * n_dot_l / (n.length() * v_light.length());
        };

        // Specular
        if s != -1 {
            let r = reflect_ray(v_light, n);
            let r_dot_v = r.dot(v);
            if r_dot_v > 0.0 {
                total_intensity += intensity * (r_dot_v / (r.length() * v.length())).powi(s);
            }
        }
        total_intensity
    })
}
