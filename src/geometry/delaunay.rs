use super::{
  circumcircle,
  io_string,
  Point,
};
use std::{
  collections::{
    BTreeMap,
    BTreeSet,
    VecDeque,
  },
  io,
};

/// Sort a, b, and c such that the resulting a, b, c represents a clockwise
/// rotation
fn sort_clockwise(points: &[Point], a: usize, b: usize, c: usize) -> (usize, usize, usize) {
  if (points[b] - points[a]).perp(&(points[c] - points[b])) < 0f32 {
    (a, b, c)
  } else {
    (a, c, b)
  }
}

/// sorts an edge a-b so that the resulting edge has the lowest-index point
/// first e.g. (5,2) becomes(2,5) but (1,3) stays the same
/// used for inserting edges into map structures (e.g. BTreeMap)
fn sort_edge_lex(a: usize, b: usize) -> (usize, usize) {
  if a > b {
    (b, a)
  } else {
    (a, b)
  }
}

/// Tests whether the circumcircle of `t` contains `point`.
fn circumcircle_contains(t: &Triad, p: Point, points: &[Point]) -> bool {
  let d = (p - t.cc_centre).norm();
  log::info!("flip: {} <= {}?", d, t.cc_radius);
  d <= t.cc_radius
}

/// returns a new triad, with 'x' as the first vertex
fn triad_normalise(triad: &Triad, x: usize) -> Triad {
  if triad.a == x {
    // nothing to do...
    triad.clone()
  } else if triad.b == x {
    Triad {
      a: triad.b,
      b: triad.c,
      c: triad.a,

      ab: triad.bc,
      bc: triad.ca,
      ca: triad.ab,
      cc_centre: triad.cc_centre,
      cc_radius: triad.cc_radius,
    }
  } else if triad.c == x {
    Triad {
      a: triad.c,
      b: triad.a,
      c: triad.b,

      ab: triad.ca,
      bc: triad.ab,
      ca: triad.bc,

      cc_centre: triad.cc_centre,
      cc_radius: triad.cc_radius,
    }
  } else {
    panic!("{} is not a vertex in {:?}", x, triad);
  }
}

/// Sets the neighbour of the edge starting at `start` in triangle `triad_id` to
/// `neighbour`.
fn bind_edge(triad_id: usize, start: usize, neighbour: usize, triads: &mut [Triad]) {
  let mut triad = &mut triads[triad_id];
  match start {
    start if start == triad.a => triad.ab = Some(neighbour),
    start if start == triad.b => triad.bc = Some(neighbour),
    start if start == triad.c => triad.ca = Some(neighbour),
    _ => panic!("{} is not a vertex of {:?}", start, triad),
  }
}

// assumes that edge ab of `first` is the shared edge,
// and that second has been normalised so that point `a` is the same in both
// triads (i.e. the shared edge is ca in the second triad)
fn do_flip(first_id: usize, second_id: usize, triads: &mut [Triad], points: &[Point]) {
  let first = triads[first_id].clone();
  let second = triads[second_id].clone();

  log::warn!(
    "flip: {} {} {} and {} {} {}",
    first.a,
    first.b,
    first.c,
    second.a,
    second.b,
    second.c,
  );

  let mut new_second = triad_normalise(&second, first.a);
  new_second.c = first.c;
  new_second.bc = Some(first_id);
  new_second.ca = first.ca;

  let (cc_centre, cc_radius) = circumcircle(new_second.a, new_second.b, new_second.c, points).unwrap();
  new_second.cc_centre = cc_centre;
  new_second.cc_radius = cc_radius;

  if let Some(neighbour) = new_second.ca {
    bind_edge(neighbour, new_second.a, second_id, triads);
  }
  triads[second_id] = new_second;

  let mut new_first = first.clone();
  new_first.a = second.b;
  new_first.ca = Some(second_id);
  new_first.ab = second.bc;

  let (cc_centre, cc_radius) = circumcircle(new_first.a, new_first.b, new_first.c, points).unwrap();
  new_first.cc_centre = cc_centre;
  new_first.cc_radius = cc_radius;

  if let Some(neighbour) = new_first.ab {
    bind_edge(neighbour, new_first.b, first_id, triads);
  }
  triads[first_id] = new_first;
}

#[derive(Clone, Debug)]
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
  cc_centre: Point,
  cc_radius: f32,
}

impl Triad {
  fn new(a: usize, b: usize, c: usize, cc_centre: Point, cc_radius: f32) -> Triad {
    Triad {
      a,
      b,
      c,
      ab: None,
      bc: None,
      ca: None,
      cc_centre,
      cc_radius,
    }
  }

  fn centre_of_mass(&self, points: &[Point]) -> Point {
    Point::new(
      (points[self.a].coords.x + points[self.b].coords.x + points[self.c].coords.x) / 3f32,
      (points[self.a].coords.y + points[self.b].coords.y + points[self.c].coords.y) / 3f32,
    )
  }
}

/// Implements the "sweephull" delaunay triangulation algorithm
///
///
/// Basic algorithm is this:
/// - Pick one of the input points as a "pivot" point and sort the points by
///   distance to the pivot
/// - Construct an initial "seed" triangle from the pivot, it's nearest
///   neighbour, and the point that forms the smallest circumcircle with those
///   two points
/// - Process the other points in order, keeping track of a convex hull of the
///   set of points that have been added. Each new point will form a triangle
///   that shares one of the edges with the convex hull Track the triangles as
///   they are created, including tracking neighbor relationships between
///   triangles
/// - Once all the points have been processed, the triangles form a valid
///   triangulation but not a Delaunay triangulation The last step is to process
///   each pair of neighbouring triangles, check if they're a valid Delaunay
///   pair (i.e., no point falls inside the circumcircle of the other three),
///   and if not, "flip" the shared edge. After this is finished, the
///   triangulation is a Delaunay triangulation.
///
/// Complexity: O(n*log(n)). One sort and two linear passes (one for the points,
/// one for the triangles)
///
/// Implementation is in three parts:
/// - Initialisation (`new`) does the sort and constructs the seed triangle
/// - A `step` function implements the processing of a single point from the
///   sorted point array. Iteration variables are stored in this struct.
/// - After the final point has been processed, `flip_step` does the final
///   flipping pass over all the triangles.
///
/// Note that pretty much all operations here use index buffers rather than
/// using the points directly This is mainly a stylistic choice on my part, but
/// it also means that there's a single point of truth for all points, which may
/// help mitigate any floating point arithmetic errors. It also means we can do
/// some operations entirely symbolically - e.g. the "normalise" operation is
/// entirely symbolic.
pub struct Delaunay {
  // iterator variables
  origin: Point,
  hull: Vec<usize>,
  sorted_points: Vec<usize>,
  next: usize,
  triads_by_edge: BTreeMap<(usize, usize), usize>,

  points: Vec<Point>,
  triads: Vec<Triad>,
}

impl Delaunay {

  // Note that this step *invalidates* some of the iterator variables,
  // because its assumed that it is called last.
  // in particular, `triads_by_edge` will be inavlidated by any flipping done in this step.
  pub fn flip_step(&mut self) {
    log::info!("{:?}", self.triads_by_edge);
    // for each adjacent pair of triangles,
    // check if they're a valid Delaunay pair, and
    // if not, flip them
    let mut visited = BTreeSet::new();
    let mut queue = VecDeque::new();
    queue.push_front(0usize);

    while let Some(tid) = queue.pop_front() {
      let triad = &self.triads[tid];
      log::info!("flip: [{}] {}, {}, {}", tid, triad.a, triad.b, triad.c,);

      // For each edge starting at vertex `start`, we may have a neighbour `next`,
      // which we need to flip with if the point `outside` is inside its circumcircle
      for (start, next, outside) in [
        (triad.a, triad.ab, triad.c),
        (triad.b, triad.bc, triad.a),
        (triad.c, triad.ca, triad.b),
      ] {
        if let Some(next) = next {
          if !visited.contains(&next) {
            log::info!("flip: discovered neighbour");
            queue.push_front(next);
            visited.insert(tid);
            continue;
          }

          // The original sweephull implementation has a separate branch here,
          // depending on which edge is shared and where that edge falls in each triangle
          // That's probably more efficient than this, but we avoid some extra code here
          // by rotating each triangle so that on the shared edge is "ab"
          // on the tid triangle and "ca" on the "next" triangle
          // (i.e., "a" is the same in both triangles)
          //
          // This means we only need to implement one "do_flip" function.
          self.triads[next] = triad_normalise(&self.triads[next], start);
          self.triads[tid] = triad_normalise(&self.triads[tid], start);
          log::info!(
            "after normalise ({}): {} {} {}",
            self.triads[next].a,
            self.triads[next].b,
            self.triads[next].c,
            start,
          );
          if circumcircle_contains(&self.triads[next], self.points[outside], &self.points) {
            do_flip(tid, next, &mut self.triads, &self.points);
            queue.push_front(tid);
            break;
          }
        }
      }
    }
  }

  pub fn step(&mut self) -> std::task::Poll<Result<(), ()>> {
    //! Attempt to add the next point to the convex hull, and add triangles to
    //! any visible edges in the hull
    log::info!("step: triads_by_edge: {:?}", self.triads_by_edge);
    if self.next >= self.sorted_points.len() {
      log::info!("flip step");
      self.flip_step();
      return std::task::Poll::Ready(Ok(()));
    }

    let pid = self.sorted_points[self.next];
    let pt = &self.points[pid];
    let hull0 = &self.points[self.hull[0]];
    let hull_pt_v = pt - hull0;

    log::info!("step {}: {}, {}", self.next, pt.coords.x, pt.coords.y);

    // visibility test: an edge (h0, h1) is visible if the angle (h0, h1, pt) is
    // counter-clockwise.
    let visible = |hid: usize, hnext: usize, hull: &[usize], points: &[Point]| -> bool {
      let h0 = &points[hull[hid]];
      let h1 = &points[hull[hnext]];
      let angle = (h1 - h0).perp(&(pt - h1));
      let visible = angle > 0f32;

      visible
    };

    let next_hid = |hid: usize, hull_len: usize| -> usize {
      if hid >= hull_len - 1 {
        0
      } else {
        hid + 1
      }
    };
    let prev_hid = |hid: usize, hull_len: usize| -> usize {
      if hid <= 0 {
        hull_len - 1
      } else {
        hid - 1
      }
    };

    log::info!("--- find first visible edge");
    // Find the origin first visible edge of the hull
    let hid_start = (0..self.hull.len())
      .rev()
      .filter(|hid| {
        let (hnext, hprev) = (
          next_hid(*hid, self.hull.len()),
          prev_hid(*hid, self.hull.len()),
        );
        visible(*hid, hnext, &self.hull, &self.points)
          && !visible(hprev, *hid, &self.hull, &self.points)
      })
      .next();

    let hid_start = match hid_start {
      None => {
        log::error!("No visible edges in the hull");
        return std::task::Poll::Ready(Err(()));
      }
      Some(v) => v,
    };

    // Any points in the hull that are not
    //  - the origin of the first visible edge (hid_start), or
    //  - the destination of the last visible edge
    // can be removed from the hull - we're about to add triangles that will make
    // them interior points.

    log::info!("--- walk visible edges");
    let mut hid = hid_start;
    loop {
      let hnext = next_hid(hid, self.hull.len());
      if !visible(hid, hnext, &self.hull, &self.points) {
        break;
      }

      log::info!("Found new visible edge");

      let h0 = self.hull[hid];
      let h1 = self.hull[hnext];

      let (centre, radius) = circumcircle(h0, h1, pid, &self.points).unwrap();
      let tid = self.triads.len();
      let triad = Triad::new(h0, h1, pid, centre, radius);
      self.triads.push(triad);
      log::info!("new triad with vertices {}, {}, {}", h0, h1, pid);

      for (a, b) in [(h0, h1), (h1, pid), (pid, h0)] {
        if let Some(other) = self.triads_by_edge.get(&sort_edge_lex(a, b)) {
          log::info!(
            "new neighbour with vertices {}, {}, {}",
            self.triads[*other].a,
            self.triads[*other].b,
            self.triads[*other].c
          );
          bind_edge(tid, a, *other, &mut self.triads);
          bind_edge(*other, b, tid, &mut self.triads);
        }
      }
      self.triads_by_edge.insert(sort_edge_lex(h0, h1), tid);
      self.triads_by_edge.insert(sort_edge_lex(h1, pid), tid);
      self.triads_by_edge.insert(sort_edge_lex(pid, h0), tid);

      if hid != hid_start {
        // previous edge must have been visible,
        // and hid is the origin of a visible edge
        // - we can remove hid from the hull
        self.hull.remove(hid);

        // note that because we're removing from hull, we don't need to increment hid.
        // we do need to check that it doesn't overflow though
        if hid == self.hull.len() {
          hid = 0
        };
      } else {
        hid = next_hid(hid, self.hull.len());
      }
    }
    self.hull.insert(hid, pid);

    log::info!(
      "New hull: \n{}",
      self
        .hull
        .iter()
        .map(|h| format!("{}, {}", self.points[*h].coords.x, self.points[*h].coords.y))
        .collect::<Vec<_>>()
        .join("\n")
    );

    self.next += 1;

    std::task::Poll::Pending
  }

  pub fn new(points: &[Point]) -> Delaunay {
    let mut sorted: Vec<usize> = (0usize..points.len()).collect();
    let p0 = sorted[0].clone();

    fn sort_radial(origin: &Point, indices: &mut [usize], points: &[Point]) {
      indices.sort_by(|p, q| {
        (points[*p] - origin)
          .norm()
          .partial_cmp(&(points[*q] - origin).norm())
          .unwrap()
      });
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

    let mut triads_by_edge = BTreeMap::new();
    triads_by_edge.insert(sort_edge_lex(a, b), 0);
    triads_by_edge.insert(sort_edge_lex(b, c), 0);
    triads_by_edge.insert(sort_edge_lex(c, a), 0);

    Delaunay {
      origin: points[0].clone(),
      hull: vec![a, b, c],
      sorted_points: sorted,
      next: 0,
      triads_by_edge,

      points: points.to_vec(),
      triads: vec![Triad::new(a, b, c, circumcircle_centre, min_radius)],
    }
  }

  pub fn gnuplot(&self, show_neighbours: bool, show_circumcircles: bool) -> io::Result<String> {
    io_string(|out| {

      writeln!(out, "set xrange [0:10]")?;
      writeln!(out, "set yrange [0:10]")?;

      writeln!(out, "$points << EOD")?;
      for (i, point) in self.points.iter().enumerate() {
        writeln!(
          out,
          "{} {} 2 {}",
          point.coords.x,
          point.coords.y,
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
          writeln!(
            out,
            "{} {} {} {}",
            f.coords.x,
            f.coords.y,
            Point::from(g - f).coords.x,
            Point::from(g - f).coords.y
          )?;
        }
      }
      writeln!(out, "EOD")?;

      if show_neighbours {
        writeln!(out, "$neighbours << EOD")?;
        for triad in &self.triads {
          for n in [triad.ab, triad.bc, triad.ca] {
            if let Some(n) = n {
              let centre = triad.centre_of_mass(&self.points);
              writeln!(
                out,
                "{} {} {} {}",
                centre.coords.x,
                centre.coords.y,
                Point::from(self.triads[n].centre_of_mass(&self.points) - centre)
                  .coords
                  .x,
                Point::from(self.triads[n].centre_of_mass(&self.points) - centre)
                  .coords
                  .y
              )?;
            }
          }
        }
        writeln!(out, "EOD")?;
      }

      if show_circumcircles {
        writeln!(out, "$circles << EOD")?;
        for triad in &self.triads {
          writeln!(out, "{} {} {}", triad.cc_centre.coords.x, triad.cc_centre.coords.y, triad.cc_radius)?;
        }
        writeln!(out, "EOD")?;
      }

      writeln!(out, "plot $points with points, \\")?;
      writeln!(out, "$edges with vectors nohead, \\")?;
      if show_neighbours {
        writeln!(out, "$neighbours with vectors, \\")?;
      }
      if show_circumcircles {
        writeln!(out, "$circles with circles, \\")?;
      }
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
