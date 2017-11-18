use std::ops::Sub;

#[derive(Debug, Clone, Copy)]
pub struct RGB {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Point { x: x, y: y, z: z }
    }
}

impl Sub for Point {
    type Output = Point;
    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

#[derive(Debug)]
pub struct Canvas {
    pub depth: i32,
    pub width: i32,
    pub height: i32,
    pub v_width: i32,
    pub v_height: i32,

    pub buffer: Vec<u8>,
}

impl Canvas {
    pub fn new(width: i32, height: i32) -> Self {
        Canvas {
            depth: 1,
            width,
            height,
            v_width: 1,
            v_height: 1,
            buffer: vec![0xff; (3 * width * height) as usize],
        }
    }

    pub fn to_viewport(&self, x: f32, y: f32) -> Point {
        Point::new(
            x * (self.v_width as f32) / (self.width as f32),
            y * (self.v_height as f32) / (self.height as f32),
            self.depth as f32,
        )
    }

    pub fn put_pixel(&mut self, x: i32, y: i32, color: RGB) {
        let c_x = self.width / 2 + x;
        let c_y = self.height / 2 - y - 1;
        if c_x < 0 || c_x >= self.width || c_y < 0 || c_y >= self.height {
            return;
        }

        let offset = (3 * (c_x + self.width * (self.width - c_y - 1))) as usize;
        self.buffer[offset] = color.b;
        self.buffer[offset + 1] = color.g;
        self.buffer[offset + 2] = color.r;
    }
}

pub fn rgb(r: u8, g: u8, b: u8) -> RGB {
    RGB { r, g, b }
}