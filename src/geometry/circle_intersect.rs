use super::{
  circumcircle,
  io_string,
  Point,
  Vector,
};
use std::{
  collections::{
    BTreeMap,
    BTreeSet,
    VecDeque,
  },
  io,
};

struct Circle {
  centre: Point,
  radius: f32,
}

impl Circle {
  fn intersection2(a: Circle, b: Circle) -> Option<(Point, Point)> {
    // First, check that the circles actually intersect
    let ab = b.centre - a.centre;
    let d = ab.norm();
    if d > (a.radius + b.radius) {
      return None;  // The circles are disjoint
    }
    if d < (a.radius - b.radius) {
      return None;  // one circle is contained within the other
    }
    if d == 0f32 && (a.radius == b.radius) {
      return None;  // the circles are coincident
    }

    // There will be two intersections, falling along a line and equidistant to it

    // distance along ab at which that line will intersect is given by
    let t = (a.radius * a.radius - b.radius * b.radius + d * d) / (2f32*d);
    // and the height from that line at which the circles will intersect:
    let u = (a.radius * a.radius - t * t).sqrt();

    let normal = Vector::new(ab.y, -ab.x);
    let normal = normal / normal.norm();

    let g = a.centre + (t / d) * ab; // translate to intersection on line ab
    
    return Some((
      g + u * normal,
      g - u * normal
    ))
  }

  fn intersection3(a: Circle, b: Circle, c: Circle) -> Option<Point> {
    let (x, y) = Circle::intersection2(a, b)?;
    if (c.centre-x).norm() == c.radius {
      Some(x)
    } else {
      Some(y)
    }
  }

  fn graphviz(&self, writer: &mut impl io::Write) {
    writeln!(writer, "{} {} {}", self.centre.coords.x, self.centre.coords.y, self.radius);
  }
}

#[test]
fn intersect2() {
  let a = Circle {
    centre: Point::new(0.0, 0.0),
    radius: 4.0,
  };
  let b = Circle {
    centre: Point::new(3.0, 0.0),
    radius: 5.0
  };

  assert_eq!(Circle::intersection2(a, b), Some((Point::new(0.0, -4.0), Point::new(0.0, 4.0))));
}

#[test]
fn intersect3() {
  let a = Circle {
    centre: Point::new(0.0, 0.0),
    radius: 4.0,
  };
  let b = Circle {
    centre: Point::new(3.0, 0.0),
    radius: 5.0
  };
  let c = Circle {
    centre: Point::new(0.0, 6.0),
    radius: 2.0,
  };

  assert_eq!(Circle::intersection3(a, b, c), Some(Point::new(0.0, 4.0)));
}