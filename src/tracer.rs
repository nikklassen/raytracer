use std::ops::Deref;
use std::f32;
use canvas::{rgb, Canvas, Point, RGB};
use rayon::prelude::*;


const INF: f32 = f32::MAX;

#[derive(Debug, Clone, Copy)]
struct ShapeProperties {
    color: RGB,
    specular: i32,
    reflective: f32,
}

trait Shape {
    fn norm(&self, p: Point, face: Object) -> Point;
    fn props(&self) -> ShapeProperties;
    fn intersect_ray(&self, origin: Point, d: Point) -> (f32, f32, Option<Object>);
}

#[derive(Debug, Clone)]
struct Sphere {
    center: Point,
    radius: f32,
    props: ShapeProperties,
}

impl Shape for Sphere {
    fn norm(&self, p: Point, _face: Object) -> Point {
        p - self.center
    }

    fn props(&self) -> ShapeProperties {
        self.props
    }

    fn intersect_ray(&self, origin: Point, d: Point) -> (f32, f32, Option<Object>) {
        let center = self.center;
        let r = self.radius;
        let oc = origin - center;

        let k1 = d.dot(d);
        let k2 = 2.0 * oc.dot(d);
        let k3 = oc.dot(oc) - r.powi(2);

        let discriminant: f32 = k2.powi(2) - 4.0 * k1 * k3;
        if discriminant < 0.0 {
            (INF, INF, None)
        } else {
            (
                ((-k2 + discriminant.sqrt()) / (2.0 * k1)),
                ((-k2 - discriminant.sqrt()) / (2.0 * k1)),
                Some(Object::Sphere(self.clone())),
            )
        }
    }
}

#[derive(Debug, Clone)]
struct Triangle {
    p1: Point,
    p2: Point,
    p3: Point,
    props: ShapeProperties,
}

/*
impl Triangle {
    fn includes_point(&self, p: Point) -> bool {
        let v0 = self.p3 - self.p1;
        let v1 = self.p2 - self.p1;
        let v2 = p - self.p1;
        let dot00 = v0.dot(v0);
        let dot01 = v0.dot(v1);
        let dot02 = v0.dot(v2);
        let dot11 = v1.dot(v1);
        let dot12 = v1.dot(v2);

        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

        u >= 0.0 && v >= 0.0 && u + v < 1.0
    }
}
*/

impl Shape for Triangle {
    fn norm(&self, p: Point, _face: Object) -> Point {
        let pa = p - self.p1;
        let pb = p - self.p2;
        let n = pb.cross(pa);
        if n.dot(p) >= 0.0 {
            n.mul(-1.0)
        } else {
            n
        }
    }

    fn props(&self) -> ShapeProperties {
        self.props
    }

    fn intersect_ray(&self, origin: Point, d: Point) -> (f32, f32, Option<Object>) {
        let e1 = self.p2 - self.p1;
        let e2 = self.p3 - self.p1;
        let h = d.cross(e2);
        let a = e1.dot(h);
        const EPS: f32 = 0.0001;
        if a > -EPS && a < EPS {
            return (INF, INF, None);
        }

        let f = 1.0/a;
        let s = origin - self.p1;
        let u = f * s.dot(h);
        if u < 0.0 || u > 1.0 {
            return (INF, INF, None);
        }

        let q = s.cross(e1);
        let v = f * d.dot(q);
        if v < 0.0 || u + v > 1.0 {
            return (INF, INF, None);
        }

        let t = f * e2.dot(q);
        if t > EPS {
            (t, t, Some(Object::Triangle(self.clone())))
        } else {
            (INF, INF, None)
        }
    }
}

#[derive(Debug, Clone)]
struct Polygon {
    faces: [Triangle; 4],
    props: ShapeProperties,
}

impl Polygon {
    pub fn new(
        p1: Point,
        p2: Point,
        p3: Point,
        p4: Point,
        props: ShapeProperties,
    ) -> Self {
        Polygon {
            faces: [
                Triangle { p1: p1, p2: p2, p3: p3, props },
                Triangle { p1: p2, p2: p3, p3: p4, props },
                Triangle { p1: p3, p2: p4, p3: p1, props },
                Triangle { p1: p4, p2: p1, p3: p2, props },
            ],
            props,
        }
    }
}

impl Shape for Polygon {
    fn norm(&self, p: Point, face: Object) -> Point {
        face.norm(p, face.clone())
    }

    fn props(&self) -> ShapeProperties {
        self.props
    }

    fn intersect_ray(&self, origin: Point, d: Point) -> (f32, f32, Option<Object>) {
        let mut closest_t = INF;
        let mut closest_face = None;
        for face in self.faces.iter() {
            let (t, _, face_obj) = face.intersect_ray(origin, d);
            if t < closest_t {
                closest_t = t;
                closest_face = face_obj;
            }
        }
        (closest_t, closest_t, closest_face)
    }
}

#[derive(Debug, Clone)]
enum Object {
    Sphere(Sphere),
    Triangle(Triangle),
    Polygon(Polygon),
}

impl Deref for Object {
    type Target = Shape;

    fn deref<'a>(&'a self) -> &'a (Shape + 'static) {
        match self {
            &Object::Sphere(ref s) => s,
            &Object::Triangle(ref t) => t,
            &Object::Polygon(ref p) => p,
        }
    }
}

#[derive(Debug, Clone)]
enum Light {
    Ambient { intensity: f32 },
    Point { intensity: f32, position: Point },
    Directional { intensity: f32, direction: Point },
}

#[derive(Debug, Clone)]
struct Scene {
    objects: Vec<Object>,
    lights: Vec<Light>,
}

const BACKGROUND_COLOR: RGB = RGB { r: 0, g: 0, b: 0 };

pub fn draw(canvas: &mut Canvas, x_rot: f32, y_rot: f32) {
    let scene: Scene = Scene {
        objects: vec![
            // Shiny red sphere
            Object::Sphere(Sphere {
                center: Point::new(0.0, -1.0, 3.0),
                radius: 1.0,
                props: ShapeProperties {
                    color: rgb(255, 0, 0),
                    specular: 500,
                    reflective: 0.0,
                },
            }),
            // Shiny blue sphere
            Object::Sphere(Sphere {
                center: Point::new(2.0, 0.0, 4.0),
                radius: 1.0,
                props: ShapeProperties {
                    color: rgb(0, 0, 255),
                    specular: 500,
                    reflective: 0.0,
                },
            }),
            // Matte green sphere
            Object::Sphere(Sphere {
                center: Point::new(-2.0, 0.0, 4.0),
                radius: 1.0,
                props: ShapeProperties {
                    color: rgb(0, 255, 0),
                    specular: 10,
                    reflective: 0.0,
                },
            }),
            // Very shiny yellow sphere
            Object::Sphere(Sphere {
                center: Point::new(0.0, -5001.0, 0.0),
                radius: 5000.0,
                props: ShapeProperties {
                    color: rgb(255, 255, 0),
                    specular: 1000,
                    reflective: 0.5,
                },
            }),
            Object::Polygon(Polygon::new(
                Point::new(-1.0, 2.0, 3.0),
                Point::new(0.0, 1.0, 3.0),
                Point::new(1.0, 2.0, 3.0),
                Point::new(0.4, 1.5, 1.0),
                ShapeProperties {
                    color: rgb(100, 100, 100),
                    specular: 10,
                    reflective: 0.0,
                },
            )),
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

    let origin = Point::new(0.0, 0.0, -1.0);
    let camera_rotation = make_rotation(x_rot / 10.0, y_rot / 10.0, 0.0);
    let segments = vec![
        [-canvas.width / 2, 0, -canvas.height / 2, 0],
        [-canvas.width / 2, 0, 0, canvas.height / 2],
        [0, canvas.width / 2, -canvas.height / 2, 0],
        [0, canvas.width / 2, 0, canvas.height / 2],
    ];

    let chunk_size = (canvas.width * canvas.height / 4) as usize;
    {
        let color_results = segments
            .par_iter()
            .map(|segment| {
                let mut result = vec![(0, 0, RGB { r: 0, g: 0, b: 0 }); chunk_size];
                let mut offset = 0;

                for x in segment[0]..segment[1] {
                    for y in segment[2]..segment[3] {
                        let p = canvas.to_viewport(x as f32, y as f32);
                        let d = point_mul(&camera_rotation, p);
                        let color = trace_ray(&scene, origin, d, canvas.depth as f32, INF, 3);
                        result[offset] = (x, y, color);
                        offset += 1;
                    }
                }
                result
            })
            .collect::<Vec<_>>();
        for chunk in color_results {
            for (x, y, color) in chunk {
                canvas.put_pixel(x, y, color);
            }
        }
    }
}

fn make_rotation(theta_x: f32, theta_y: f32, theta_z: f32) -> Vec<Vec<f32>> {
    let rx = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, theta_x.cos(), -theta_x.sin()],
        vec![0.0, theta_x.sin(), theta_x.cos()],
    ];
    let ry = vec![
        vec![theta_y.cos(), 0.0, theta_y.sin()],
        vec![0.0, 1.0, 0.0],
        vec![-theta_y.sin(), 0.0, theta_y.cos()],
    ];
    let rz = vec![
        vec![theta_z.cos(), -theta_z.sin(), 0.0],
        vec![theta_z.sin(), theta_z.cos(), 0.0],
        vec![0.0, 0.0, 1.0],
    ];
    let m = matrix_mul2(&ry, &rx);
    matrix_mul2(&rz, &m)
}

fn matrix_mul2(m1: &Vec<Vec<f32>>, m2: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let n = m1.len();
    let m = m2[0].len();
    let mut result: Vec<Vec<f32>> = vec![vec![0.0; m]; n];
    for row in 0..n {
        for col in 0..m {
            for i in 0..n {
                result[row][col] += m1[row][i] * m2[i][col];
            }
        }
    }
    result
}

fn point_mul(m1: &Vec<Vec<f32>>, p: Point) -> Point {
    Point {
        x: m1[0][0] * p.x + m1[0][1] * p.y + m1[0][2] * p.z,
        y: m1[1][0] * p.x + m1[1][1] * p.y + m1[1][2] * p.z,
        z: m1[2][0] * p.x + m1[2][1] * p.y + m1[2][2] * p.z,
    }
}

fn closest_intersection<'a>(
    scene: &Scene,
    origin: Point,
    d: Point,
    t_min: f32,
    t_max: f32,
) -> (Option<Object>, f32, Option<Object>) {
    let mut closest_t = INF;
    let mut closest_sphere = None;
    let mut closest_face = None;
    for obj in scene.objects.iter() {
        let (t1, t2, face) = obj.intersect_ray(origin, d);
        if t1 >= t_min && t1 <= t_max && t1 < closest_t {
            closest_t = t1;
            closest_sphere = Some(obj.clone());
            closest_face = face.clone();
        }
        if t2 >= t_min && t2 <= t_max && t2 < closest_t {
            closest_t = t2;
            closest_sphere = Some(obj.clone());
            closest_face = face.clone();
        }
    }
    (closest_sphere, closest_t, closest_face)
}

fn has_shadow(scene: &Scene, origin: Point, d: Point, t_min: f32, t_max: f32) -> bool {
    // similar to find closest intersection, except we just find if there is an intersection
    for obj in scene.objects.iter() {
        let (t1, t2, _) = obj.intersect_ray(origin, d);
        if t1 >= t_min && t1 <= t_max && t1 < INF || t2 >= t_min && t2 <= t_max && t2 < INF {
            return true;
        }
    }
    false
}

fn trace_ray(scene: &Scene, origin: Point, d: Point, t_min: f32, t_max: f32, depth: u32) -> RGB {
    let (closest_sphere, closest_t, face) = closest_intersection(scene, origin, d, t_min, t_max);
    if closest_sphere.is_none() {
        return BACKGROUND_COLOR;
    }
    let sphere = closest_sphere.unwrap();

    // Compute local color
    let p = origin + d.mul(closest_t);
    let mut n = sphere.norm(p, face.unwrap());
    n = n.div(n.length());
    let props = sphere.props();
    let local_color = props
        .color
        .with_intensity(compute_lighting(scene, p, n, d.mul(-1.0), props.specular));

    // If we hit the recursion limit or the object is not relective, we're done
    let r = props.reflective;
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
        if has_shadow(scene, p, v_light, 0.001, t_max) {
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
