#![allow(dead_code, unused_variables)]

use nalgebra::{
  self as na,
  Point2,
  Vector2,
};
use std::{
  cmp::Ordering,
  io,
};

mod delaunay;
pub use delaunay::*;

#[cfg(feature = "voronoi")]
mod voronoi;
#[cfg(feature = "voronoi")]
pub use voronoi::*;

pub type Point = Point2<f32>;

pub fn point_cmp_lex_xy(lhs: &Point, rhs: &Point) -> Option<Ordering> {
  match lhs.coords.x.partial_cmp(&rhs.coords.x) {
    Some(Ordering::Equal) => lhs.coords.y.partial_cmp(&rhs.coords.y),
    other => other,
  }
}

pub fn point_cmp_lex_yx(lhs: &Point, rhs: &Point) -> Option<Ordering> {
  match lhs.coords.y.partial_cmp(&rhs.coords.y) {
    Some(Ordering::Equal) => lhs.coords.x.partial_cmp(&rhs.coords.x),
    other => other,
  }
}

fn circumcircle(a: usize, b: usize, c: usize, points: &[Point]) -> Option<(Point, f32)> {
  //! Returns the centrepoint and radius of the circumcircle of three points
  use nalgebra::Matrix2;

  // We can get the centrepoint by finding the:
  //  - edges of the triangle
  //  - normals at their midpoints
  //  - their intersection
  let (a, b, c) = (points[a], points[b], points[c]);
  // midpoints
  let (mab, mbc, mca) = (na::center(&b, &a), na::center(&c, &b), na::center(&a, &c));
  // normals
  let (nab, nbc, nca) = (
    Vector2::new((a - b).y, (b - a).x),
    Vector2::new((b - c).y, (c - b).x),
    Vector2::new((c - a).y, (a - c).x),
  );

  let intersection = if let Some(mat) = Matrix2::from_columns(&[nab, nbc]).try_inverse() {
    // The column matrix having an inverse means the simultaneous equations
    // for the two normals can be solved.
    // Their solution in terms of scales of nab and nbc is given by
    let sols = Point2::from(mat * (mbc - mab));
    mab + nab * sols.coords[0]
  } else if let Some(mat) = Matrix2::from_columns(&[nbc, nca]).try_inverse() {
    let sols = Point2::from(mat * (mca - mbc));
    mbc + nbc * sols.coords[0]
  } else if let Some(mat) = Matrix2::from_columns(&[nca, nab]).try_inverse() {
    let sols = Point2::from(mat * (mab - mca));
    mca + nca * sols.coords[0]
  } else {
    log::warn!("points {}, {}, and {} are colinear", mab, mbc, mca);
    return None;
  };

  let radius = (intersection - a).norm();
  Some((intersection, radius))
}

#[cfg(test)]
#[test]
fn circumcircles() {
  let a = Point::new(0f32, 1f32);
  let b = Point::new(1f32, 2f32);
  let c = Point::new(2f32, 1f32);
  assert_eq!(
    circumcircle(0, 1, 2, &[a, b, c]),
    Some((Point::new(1f32, 1f32), 1f32))
  );
}

#[derive(PartialEq, Debug)]
pub struct Box2 {
  pub tl: Point,
  pub tr: Point,
  pub br: Point,
  pub bl: Point,
}

impl Box2 {
  pub fn as_points(&self) -> Vec<Point> { vec![self.tl, self.tr, self.br, self.bl] }
}

pub struct Bounds {
  pub top: f32,
  pub bottom: f32,
  pub left: f32,
  pub right: f32,
}

impl Into<Box2> for Bounds {
  fn into(self) -> Box2 {
    Box2 {
      tl: Point::new(self.left, self.top),
      tr: Point::new(self.right, self.top),
      br: Point::new(self.right, self.bottom),
      bl: Point::new(self.left, self.bottom),
    }
  }
}

impl Box2 {
  pub fn bounding_box(points: &[Point]) -> Option<Box2> {
    if points.is_empty() {
      return None;
    }
    let mut bounds = Bounds {
      top: points[0].coords.y,
      bottom: points[0].coords.y,
      left: points[0].coords.x,
      right: points[0].coords.x,
    };
    for point in points {
      bounds = Bounds {
        top: bounds.top.min(point.y),
        bottom: bounds.bottom.max(point.y),
        left: bounds.left.min(point.x),
        right: bounds.right.max(point.x),
      }
    }

    Some(bounds.into())
  }
}

/// Computes the convex hull of a set of points by the monotone chain
/// algorithm see https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
///
/// Algorithm: Assume the list is sorted by x, then y
/// Build the "upper" and "lower" sections of the hull separately,
/// where "upper" and "lower" mean above and below the vertical position of
/// the first point.
///
/// To build each hull, keep a stack of points currently considered part of
/// the hull Attempt to join each new point to the hull by considering the
/// preceeding two points in the hull and the new point If the angle is
/// concave, great! add the new point If it's convex, drop the current last
/// point in the hull and retry If there are less than two points in the
/// hull (so no angle can be calculated), just add the new point.
///
/// To compute the lower hull we do the same as the upper, but traverse the
/// points in reverse. This has the effect of inverting the
/// perpendicular_product computation.
///
/// Returns a sequence of indices into `points` giving the order of points
/// included in the convex hull
pub fn convex_hull(points: &[Point]) -> Vec<usize> {
  fn perpendicular_product(a: Point, o: Point, b: Point) -> f32 { (a - o).perp(&(b - o)) }

  fn compute_hull<'a>(it: impl Iterator<Item = &'a usize>, points: &[Point]) -> Vec<usize> {
    let mut hull = Vec::new();
    for pi in it {
      while hull.len() >= 2
        && perpendicular_product(
          points[hull[hull.len() - 2]],
          points[hull[hull.len() - 1]],
          points[*pi],
        ) <= 0f32
      {
        hull.pop();
      }
      hull.push(*pi);
    }

    // The last element of each hull is the first element of the next,
    // so we drop it here to avoid duplicate points in the merged hull
    hull.pop();

    return hull;
  }

  let mut sorted = (0..points.len()).collect::<Vec<usize>>();
  sorted.sort_by(|a, b| point_cmp_lex_xy(&points[*a], &points[*b]).unwrap());

  compute_hull(sorted.iter(), points)
    .into_iter()
    .chain(compute_hull(sorted.iter().rev(), points))
    .collect::<Vec<_>>()
}

fn io_string(mut f: impl FnMut(&mut dyn io::Write) -> io::Result<()>) -> io::Result<String> {
  use io::{
    Read,
    Seek,
  };

  let mut out = io::Cursor::new(Vec::new());

  f(&mut out)?;

  let mut result = String::new();
  out.seek(io::SeekFrom::Start(0))?;
  out.read_to_string(&mut result)?;
  Ok(result)
}

#[cfg(test)]
use rand::{
  prelude::SliceRandom,
  thread_rng,
};

#[cfg(test)]
#[test]
fn bounding_boxes() {
  let mut points = Vec::new();

  let bounds = Box2 {
    tl: Point::new(0f32, 0f32),
    tr: Point::new(10f32, 0f32),
    br: Point::new(10f32, 10f32),
    bl: Point::new(0f32, 10f32),
  };

  for x in 0..=10 {
    for y in 0..=10 {
      points.push(Point::new(x as f32, y as f32));
    }
  }

  for _ in 0..100 {
    points.shuffle(&mut thread_rng());
    assert_eq!(Box2::bounding_box(&points).unwrap(), bounds);
  }
}

#[cfg(test)]
#[test]
fn convex_hulls_fixed() {
  let mut points = Vec::new();

  let expected = vec![
    Point::new(0f32, 0f32),
    Point::new(0f32, 2f32),
    Point::new(2f32, 2f32),
    Point::new(2f32, 0f32),
  ];

  for x in 0..=2 {
    for y in 0..=2 {
      points.push(Point::new(x as f32, y as f32));
    }
  }

  for _ in 0..4 {
    points.shuffle(&mut thread_rng());
    assert_eq!(
      convex_hull(&points)
        .iter()
        .map(|i| points[*i])
        .collect::<Vec<_>>(),
      expected
    );
  }
}
