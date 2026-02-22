use glam::{DVec2, IVec2};
use serde::{Deserialize, Serialize};

use crate::{
    Color,
    img::{Blending, PixelData, RawImage},
    librt2d::{Hit2d, Object, Ray2d},
};

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct Aabb {
    pub center: DVec2,
    pub half_size: DVec2,
}
impl Aabb {
    pub fn new(center: DVec2, half_size: DVec2) -> Self {
        Self { center, half_size }
    }

    pub fn contains(&self, other: &Self) -> bool {
        (self.center - self.half_size)
            .cmplt(other.center - other.half_size)
            .all()
            && (self.center + self.half_size)
                .cmpgt(other.center + other.half_size)
                .all()
    }

    pub fn union(&self, other: &Self) -> Self {
        let min = (self.center - self.half_size).min(other.center - other.half_size);
        let max = (self.center + self.half_size).max(other.center + other.half_size);

        let center = min.midpoint(max);
        let half_size = (max - min) / 2.;

        Aabb::new(center, half_size)
    }

    pub fn hit(&self, r: &Ray2d) -> Option<f64> {
        // https://tavianator.com/2011/ray_box.html
        let tx1 = (self.center.x - self.half_size.x - r.origin.x) * r.dir.recip().x;
        let tx2 = (self.center.x + self.half_size.x - r.origin.x) * r.dir.recip().x;

        let tmin = tx1.min(tx2);
        let tmax = tx1.max(tx2);

        let ty1 = (self.center.y - self.half_size.y - r.origin.y) * r.dir.recip().y;
        let ty2 = (self.center.y + self.half_size.y - r.origin.y) * r.dir.recip().y;

        let tmin = ty1.min(ty2).max(tmin);
        let tmax = ty1.max(ty2).min(tmax);

        //if tmax > tmin && tmax > 0. {
        //    Some(tmax)
        //} else {
        //    None
        //}
        if tmax >= tmin.max(0.0) {
            Some(tmin.max(0.0))
        } else {
            None
        }
    }

    pub fn draw(&self, raw_image: &mut RawImage, color: Color) {
        let blend = Blending::Replace;

        //// cross
        //for i in -3..=3 {
        //    let pixel = self.center.as_ivec2() + IVec2::new(i, 0);
        //    let _ = raw_image.draw_pixel(pixel, color, blend);
        //}
        //for j in -3..=3 {
        //    let pixel = self.center.as_ivec2() + IVec2::new(0, j);
        //    let _ = raw_image.draw_pixel(pixel, color, blend);
        //}

        let pixel_data = PixelData {
            value: color,
            weight: 1.,
        };

        // box TOP
        let h = self.half_size.x as i32;
        for i in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(i, -self.half_size.y as i32);
            let _ = raw_image.draw_pixel(pixel, pixel_data, blend);
        }
        // box BOTTOM
        for i in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(i, self.half_size.y as i32);
            let _ = raw_image.draw_pixel(pixel, pixel_data, blend);
        }
        // box LEFT
        let h = self.half_size.y as i32;
        for j in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(self.half_size.x as i32, j);
            let _ = raw_image.draw_pixel(pixel, pixel_data, blend);
        }
        // box RIGHT
        for j in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(-self.half_size.x as i32, j);
            let _ = raw_image.draw_pixel(pixel, pixel_data, blend);
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Node {
    objects: Vec<Object>,
    pub aabb: Aabb,
    pub top_right: Option<Box<Node>>,
    pub top_left: Option<Box<Node>>,
    pub bottom_left: Option<Box<Node>>,
    pub bottom_right: Option<Box<Node>>,
}

#[derive(Copy, Clone)]
pub enum Direction {
    TopRight,
    TopLeft,
    BottomLeft,
    BottomRight,
}

impl Node {
    pub fn new(aabb: Aabb) -> Self {
        let objects = Vec::new();
        Self {
            objects,
            aabb,
            top_right: None,
            top_left: None,
            bottom_left: None,
            bottom_right: None,
        }
    }

    pub fn insert(&mut self, obj: Object) -> bool {
        use Direction::*;
        let obj_aabb = obj.aabb();

        // too big to fit in children
        if obj_aabb.half_size.cmpgt(self.aabb.half_size / 2.).any() {
            self.objects.push(obj);
            return true;
        }

        let sub_aabbs = [BottomLeft, BottomRight, TopLeft, TopRight].map(|dir| self.sub_aabb(dir));

        for (node, sub_aabb) in self.iter_mut().iter_mut().zip(sub_aabbs) {
            if sub_aabb.contains(&obj_aabb) {
                let child = node.get_or_insert_with(|| Box::new(Node::new(sub_aabb)));
                return child.insert(obj);
            }
        }
        // does not fit in any children

        if self.aabb.contains(&obj_aabb) {
            self.objects.push(obj);
            return true;
        }

        // can't insert
        return false;
    }

    pub fn sub_aabb(&self, dir: Direction) -> Aabb {
        let half_size = self.aabb.half_size / 2.;
        let center = match dir {
            Direction::TopRight => self.aabb.center + DVec2::new(half_size.x, half_size.y),
            Direction::TopLeft => self.aabb.center + DVec2::new(-half_size.x, half_size.y),
            Direction::BottomLeft => self.aabb.center + DVec2::new(-half_size.x, -half_size.y),
            Direction::BottomRight => self.aabb.center + DVec2::new(half_size.x, -half_size.y),
        };
        Aabb::new(center, half_size)
    }

    pub fn iter(&self) -> [&Option<Box<Node>>; 4] {
        [
            &self.bottom_left,
            &self.top_left,
            &self.top_right,
            &self.bottom_right,
        ]
    }

    pub fn iter_mut(&mut self) -> [&mut Option<Box<Node>>; 4] {
        [
            &mut self.bottom_left,
            &mut self.top_left,
            &mut self.top_right,
            &mut self.bottom_right,
        ]
    }

    pub fn hit(&self, ray: &Ray2d) -> Option<(&Object, Hit2d)> {
        if self.aabb.hit(&ray).is_none() {
            // ray does not hit current node nor children
            return None;
        }

        // check if ray hits children
        let children_hits: Vec<(&Object, Hit2d)> = [
            &self.bottom_left,
            &self.bottom_right,
            &self.top_left,
            &self.top_right,
        ]
        .iter()
        .filter_map(|node| {
            node.as_ref().map(|node| {
                //
                match node.aabb.hit(ray) {
                    None => None,
                    Some(_) => node.hit(ray), // recurse
                }
            })
        })
        .flatten()
        .collect();

        // object could be in current node?
        let local_hits: Vec<(&Object, Hit2d)> = self
            .objects
            .iter()
            .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj, hit)))
            .collect();

        let mut hits = Vec::new();
        hits.extend(children_hits);
        hits.extend(local_hits);

        hits.sort_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));

        hits.first().cloned()
    }

    pub fn draw(&self, raw_image: &mut RawImage, color: Color) {
        self.aabb.draw(raw_image, color);

        for sub_node in [
            &self.top_right,
            &self.top_left,
            &self.bottom_right,
            &self.bottom_left,
        ] {
            match sub_node {
                Some(sub_node) => sub_node.draw(raw_image, color),
                None => (),
            }
        }
    }
}
