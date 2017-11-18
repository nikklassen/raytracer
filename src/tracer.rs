use std::f32;
use canvas::{RGB, Point, Canvas, rgb};

const INF: f32 = f32::MAX;

#[derive(Debug, Clone)]
struct Sphere {
    center: Point,
    radius: f32,
    color: RGB,
}

#[derive(Debug)]
struct Scene {
    spheres: Vec<Sphere>,
}

const BACKGROUND_COLOR: RGB = RGB {
    r: 255,
    g: 255,
    b: 255,
};


pub fn draw(canvas: &mut Canvas) {
    let origin = Point::new(0.0, 0.0, 0.0);
    for x in (-canvas.width / 2)..(canvas.width / 2) {
        for y in (-canvas.height / 2)..(canvas.height / 2) {
            let d = canvas.to_viewport(x as f32, y as f32);
            let color = trace_ray(origin, d, canvas.depth as f32, INF);
            canvas.put_pixel(x, y, color);
        }
    }
}

fn trace_ray(origin: Point, d: Point, t_min: f32, t_max: f32) -> RGB {
    let mut closest_t = INF;
    let mut closest_sphere = None;
    let scene = Scene {
        spheres: vec![
            Sphere {
                center: Point::new(0.0, -1.0, 3.0),
                radius: 1.0,
                color: rgb(255, 0, 0),
            },
            Sphere {
                center: Point::new(2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 0, 255),
            },
            Sphere {
                center: Point::new(-2.0, 0.0, 4.0),
                radius: 1.0,
                color: rgb(0, 255, 0),
            },
        ],
    };
    for sphere in scene.spheres {
        let (t1, t2) = intersect_ray_sphere(origin, d, &sphere);
        if t1 >= t_min && t1 <= t_max && t1 < closest_t {
            closest_t = t1;
            closest_sphere = Some(sphere.clone());
        }
        if t2 >= t_min && t2 <= t_max && t2 < closest_t {
            closest_t = t2;
            closest_sphere = Some(sphere.clone());
        }
    }
    if closest_sphere.is_none() {
        BACKGROUND_COLOR
    } else {
        closest_sphere.unwrap().color
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