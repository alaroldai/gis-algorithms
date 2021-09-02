use {
  super::{Point, circumcircle, io_string},
  std::{
    io,
    collections::VecDeque,
  },
};

// Sort a, b, and c such that the resulting a, b, c represents a clockwise rotation
fn sort_clockwise(points: &[Point], a: usize, b: usize, c: usize) -> (usize, usize, usize) {
  if (points[b]-points[a]).perp(&(points[c]-points[b])) < 0f32 {
    (a, b, c)
  } else {
    (a, c, b)
  }
}

struct HullPoint {
  point: usize,
  triad: usize,
}

pub struct Triad {
  // points
  a: usize,
  b: usize,
  c: usize,

  // triangle ids for neighbors
  ab: Option<usize>,
  bc: Option<usize>,
  ca: Option<usize>,

  // circumcircle
  centre: Point,
  radius: f32,
}

impl Triad {
  fn new(a: usize, b: usize, c: usize, centre: Point, radius: f32) -> Triad {
    Triad {
      a, b, c,
      ab: None,
      bc: None,
      ca: None,
      centre, radius,
    }
  }
}

pub struct Delaunay {
  // iterator variables
  origin: Point,
  hull: Vec<usize>,
  sorted_points: Vec<usize>,
  next: usize,

  points: Vec<Point>,
  triads: Vec<Triad>,
}

impl Delaunay {

  pub fn step(&mut self) -> std::task::Poll<Result<(), ()>> {
    //! Attempt to add the next point to the convex hull, and add triangles to any visible edges in the hull
    if self.next >= self.sorted_points.len() {
      return std::task::Poll::Ready(Ok(()));
    }

    let pid = self.sorted_points[self.next];
    let pt = &self.points[pid];
    let hull0 = &self.points[self.hull[0]];
    let hull_pt_v = pt - hull0;

    log::info!("step {}: {}, {}", self.next, pt.coords.x, pt.coords.y);

    log::info!("Hull: \n{}", self.hull.iter().map(|h| format!("{}, {}", self.points[*h].coords.x, self.points[*h].coords.y)).collect::<Vec<_>>().join("\n"));

    // visibility test: an edge (h0, h1) is visible if the angle (h0, h1, pt) is counter-clockwise.
    let visible = |hid: usize, hnext: usize, hull: &[usize], points: &[Point]| -> bool {
      let h0 = &points[hull[hid]];
      let h1 = &points[hull[hnext]];
      let angle = (h1-h0).perp(&(pt-h1));
      let visible = angle > 0f32;
      log::info!(
        "Edge {}: ({}, {}) -> {}: ({}, {}) is {} ",
        hid, h0.coords.x, h0.coords.y,
        hnext, h1.coords.x, h1.coords.y,
        if visible { "visible" } else { "not visible" }
      );
      visible
    };

    let next_hid = |hid: usize, hull_len: usize| -> usize { if hid >= hull_len - 1 { 0 } else { hid + 1} };
    let prev_hid = |hid: usize, hull_len: usize| -> usize { if hid <= 0 { hull_len - 1 } else { hid - 1} };

    log::info!("--- find first visible edge");
    // Find the origin first visible edge of the hull
    let hid_start = (0..self.hull.len()).rev().filter(|hid| {
      let (hnext, hprev) = (next_hid(*hid, self.hull.len()), prev_hid(*hid, self.hull.len()));
      visible(*hid, hnext, &self.hull, &self.points)
      && !visible(hprev, *hid, &self.hull, &self.points)
    }).next();

    let hid_start = match hid_start {
      None => {
        log::error!("No visible edges in the hull");
        return std::task::Poll::Ready(Err(()));
      }
      Some(v) => v
    };

    // Any points in the hull that are not
    //  - the origin of the first visible edge (hid_start), or
    //  - the destination of the last visible edge
    // can be removed from the hull - we're about to add triangles that will make them interior points.

    log::info!("--- walk visible edges");
    let mut hid = hid_start;
    loop {
      let hnext = next_hid(hid, self.hull.len());
      if !visible(hid, hnext, &self.hull, &self.points) { break; }

      let h0 = self.hull[hid];
      let h1 = self.hull[hnext];
      
      let (centre, radius) = circumcircle(h0, h1, pid, &self.points).unwrap();
      self.triads.push(Triad::new(h0, h1, pid, centre, radius));
      if hid != hid_start {
        // previous edge must have been visible,
        // and hid is the origin of a visible edge
        // - we can remove hid from the hull
        log::info!("remove edge {} from hull", hid);
        self.hull.remove(hid);

        // note that because we're removing from hull, we don't need to increment hid.
        // we do need to check that it doesn't overflow though
        if hid == self.hull.len() { hid = 0 };
      } else {
        hid = next_hid(hid, self.hull.len());
      }
    }
    self.hull.insert(hid, pid);

    log::info!("New hull: \n{}", self.hull.iter().map(|h| format!("{}, {}", self.points[*h].coords.x, self.points[*h].coords.y)).collect::<Vec<_>>().join("\n"));

    self.next += 1;

    std::task::Poll::Pending
  }

  pub fn new(points: &[Point]) -> Delaunay {
    let mut sorted: Vec<usize> = (0usize..points.len()).collect();
    let p0 = sorted[0].clone();

    fn sort_radial(origin: &Point, indices: &mut [usize], points: &[Point]) {
      indices.sort_by(|p, q| (points[*p]-origin).norm().partial_cmp(&(points[*q]-origin).norm()).unwrap());
    }
    sort_radial(&points[sorted[0]], &mut sorted, points);

    // find point to complete the smallest circumcircle with 0 and 1
    let mut min_radius = f32::MAX;
    let mut circumcircle_centre = None;
    let mut c_idx = None; // index into sorted
    for i in 2..sorted.len() {
      if let Some((centre, radius)) = circumcircle(sorted[0], sorted[1], sorted[i], points) {
        if radius < min_radius {
          circumcircle_centre = Some(centre);
          min_radius = radius;
          c_idx = Some(i);
        }
      }
    }
    let (circumcircle_centre, c_idx) = (circumcircle_centre.unwrap(), c_idx.unwrap());

    // Ensure hull is sorted in clockwise order
    let (a, b, c) = sort_clockwise(&points, sorted[0], sorted[1], sorted[c_idx]);

    sorted.remove(c_idx);
    sorted.remove(1);
    sorted.remove(0);

    Delaunay {
      origin: points[0].clone(),
      hull: vec![ a, b, c, ],
      sorted_points: sorted,
      next: 0,

      points: points.to_vec(),
      triads: vec![Triad::new(a, b, c, circumcircle_centre, min_radius)],
    }
  }

  pub fn gnuplot(&self) -> io::Result<String> {
    use io::{ Write, };
    io_string(|out| {
      writeln!(out, "$points << EOD")?;
      for (i, point) in self.points.iter().enumerate() {
        writeln!(out, "{} {} 2 {}", point.coords.x, point.coords.y,
          if i == 0 { "2" } else { "6" }
        )?;
      }
      writeln!(out, "EOD")?;
      writeln!(out, "$edges << EOD")?;
      for triad in &self.triads {
        let a = self.points[triad.a];
        let b = self.points[triad.b];
        let c = self.points[triad.c];

        for (f, g) in [(a, b), (b, c), (c, a)] {
          writeln!(out, "{} {} {} {}",
            f.coords.x, f.coords.y,
            Point::from(g-f).coords.x, Point::from(g-f).coords.y
          )?;
        }
      }
      writeln!(out, "EOD")?;
      writeln!(out, "plot $points with points, \\")?;
      writeln!(out, "$edges with vectors nohead, \\")?;
      writeln!(out)?;

      Ok(())
    })
  }
}

#[cfg(test)]
#[test]
fn point_sort_clockwise() {
  let points = [
    Point::new(2f32, 1f32),
    Point::new(1f32, 2f32),
    Point::new(0f32, 0f32),
  ];

  assert_eq!(sort_clockwise(&points, 0, 1, 2), (0, 2, 1));
  assert_eq!(sort_clockwise(&points, 0, 2, 1), (0, 2, 1));

  let points = [
    Point::new(8.608241, 7.5718346),
    Point::new(9.205788, 8.384809),
    Point::new(8.754458, 9.243713),
  ];

  assert_eq!(sort_clockwise(&points, 0, 1, 2), (0, 2, 1));
  assert_eq!(sort_clockwise(&points, 0, 2, 1), (0, 2, 1));
}