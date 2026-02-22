use macroquad::prelude::*;
use raytracer::{
    aabb::{Aabb, Node},
    librt2d::{Circle, Object, Ray2d, Shape},
    spectrum::Spectrum,
};

fn draw_aabb(node: &Node, request: &Ray2d) {
    let aabb = node.aabb;
    let top_left = aabb.center - aabb.half_size * 0.99;

    let color = if node.hit(request).is_some() {
        RED
    } else {
        ORANGE
    };

    draw_rectangle_lines(
        top_left.x as f32,
        top_left.y as f32,
        aabb.half_size.x as f32 * 2. * 0.99,
        aabb.half_size.y as f32 * 2. * 0.99,
        1.,
        color,
    );

    for sub in node.iter() {
        match sub {
            None => continue,
            Some(sub_node) => draw_aabb(sub_node, request),
        }
    }
}

#[macroquad::main("BasicShapes")]
async fn main() {
    let mut quadtree = Node::new(Aabb::new(
        DVec2::new(800. / 2., 600. / 2.),
        DVec2::new(800. / 2., 600. / 2.),
    ));

    let mut objects = Vec::new();
    let mut request = Ray2d::new(DVec2::ZERO, DVec2::ZERO);
    let mut i = 0;

    loop {
        if is_mouse_button_pressed(MouseButton::Left) {
            let p = mouse_position();
            let center = DVec2::new(p.0 as f64, p.1 as f64);
            let r = 2.;
            let circle = Shape::Circle(Circle::new(center, r));
            let obj = Object::new(
                circle,
                raytracer::librt2d::Material::Emissive {
                    emission: Spectrum::default(),
                },
            );
            objects.push(obj.clone());

            let res = quadtree.insert(obj);
        }

        if is_mouse_button_pressed(MouseButton::Right) {
            let p = mouse_position();
            let p = DVec2::new(p.0 as f64, p.1 as f64);

            if i % 2 == 0 {
                request.origin = p;
            } else {
                request.dir = (p - request.origin).normalize_or_zero();
            }
            i += 1;
        }

        clear_background(BLACK);

        for obj in &objects {
            match obj.shape {
                Shape::Circle(c) => {
                    draw_circle_lines(c.center.x as f32, c.center.y as f32, c.r as f32, 1., RED)
                }
                _ => (),
            }
        }
        draw_aabb(&quadtree, &request);

        let p = request.origin;
        let p2 = request.origin + request.dir * 200.;
        draw_circle(p.x as f32, p.y as f32, 5., BLUE);
        draw_line(p.x as f32, p.y as f32, p2.x as f32, p2.y as f32, 1., BLUE);

        next_frame().await
    }
}
