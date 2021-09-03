use {
  super::{Point, circumcircle, io_string},
  std::{
    io,
    collections::{VecDeque, BTreeMap, BTreeSet},
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

fn sort_edge_lex(a: usize, b: usize) -> (usize, usize) {
  if a > b { (b, a) }
  else { (a, b) }
}

fn circumcircle_contains(t: &Triad, p: Point, points: &[Point]) -> bool {
  log::info!("testing neighbour: {}, {}, {}", t.a, t.b, t.c);
  if let Some((centre, radius)) = circumcircle(t.a, t.b, t.c, points) {
    let d = (p-centre).norm();
    log::info!("flip: {} <= {}?", d, radius);
    d <= radius
  } else {
    false
  }
}

// returns a new triad, with 'x' as the first vertex
fn triad_normalise(triad: &Triad, x: usize) -> Triad {
  if triad.a == x {
    // nothing to do...
    triad.clone()
  }
  else if triad.b == x {
    Triad {
      a: triad.b, 
      b: triad.c, 
      c: triad.a,

      ab: triad.bc,
      bc: triad.ca,
      ca: triad.ab,
      centre: triad.centre,
      radius: triad.radius,
    }
  }
  else if triad.c == x {
    Triad {
      a: triad.c, 
      b: triad.a, 
      c: triad.b,

      ab: triad.ca,
      bc: triad.ab,
      ca: triad.bc,

      centre: triad.centre,
      radius: triad.radius,
    }
  } else {
    panic!("{} is not a vertex in {:?}", x, triad);
  }
}

fn bind_edge(triad_id: usize, start: usize, neighbour: usize, triads: &mut [Triad]) {
  let mut triad = &mut triads[triad_id];
  if triad.a == start {
    triad.ab = Some(neighbour);
    log::info!("binding {}.{} = {}", triad_id, "ab", neighbour);
  }
  else if triad.b == start {
    triad.bc = Some(neighbour);
    log::info!("binding {}.{} = {}", triad_id, "bc", neighbour);
  }
  else if triad.c == start {
    triad.ca = Some(neighbour);
    log::info!("binding {}.{} = {}", triad_id, "ca", neighbour);
  }
  else { panic!("{} is not a vertex of {:?}", start, triad) }
}

// assumes that edge ab of `first` is the shared edge,
// and that second has been normalised so that point `a` is the same in both triads
// (i.e. the shared edge is ca in the second triad)
fn do_flip(first_id: usize, second_id: usize, triads: &mut [Triad], points: &[Point]) {
  let first = triads[first_id].clone();
  let second = triads[second_id].clone();

  log::info!("flip: {} {} {} and {} {} {}",
    first.a, first.b, first.c,
    second.a, second.b, second.c,
  );

  let mut new_second = triad_normalise(&second, first.a);
  new_second.c = first.c;
  new_second.bc = Some(first_id);
  new_second.ca = first.ca;

  if let Some(neighbour) = new_second.ca {
    bind_edge(neighbour, new_second.a, second_id, triads);
  }
  triads[second_id] = new_second;

  let mut new_first = first.clone();
  new_first.a = second.b;
  new_first.ca = Some(second_id);
  new_first.ab = second.bc;
  if let Some(neighbour) = new_first.ab {
    bind_edge(neighbour, new_first.b, first_id, triads);
  }
  triads[first_id] = new_first;
}

struct HullPoint {
  point: usize,
  triad: usize,
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

  fn com(&self, points: &[Point]) -> Point {
    Point::new(
      (points[self.a].coords.x + points[self.b].coords.x + points[self.c].coords.x) / 3f32,
      (points[self.a].coords.y + points[self.b].coords.y + points[self.c].coords.y) / 3f32,
    )
  }
}

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

  pub fn flip_step(&mut self) {
    log::info!("{:?}", self.triads_by_edge);
    // for each adjacent pair of triangles,
    // check if they're a valid Delaunay pair, and
    // if not, flip them
    let mut visited = BTreeSet::new();
    let mut queue = VecDeque::new();
    queue.push_front(0usize);
    
    while let Some(tid) = queue.pop_front() {
      visited.insert(tid);
      let triad = &self.triads[tid];
      log::info!("flip: [{}] {}, {}, {}", tid, triad.a, triad.b, triad.c, );

      // For each edge starting at vertex `start`, we may have a neighbour `next`,
      // which we need to flip with if the point `outside` is inside its circumcircle
      for (start, next, outside) in [
        (triad.a, triad.ab, triad.c),
        (triad.b, triad.bc, triad.a),
        (triad.c, triad.ca, triad.b)]
      {
        if let Some(next) = next {
          if !visited.contains(&next) {
            log::info!("flip: discovered neighbour");
            queue.push_back(next);
            continue;
          }

          self.triads[next] = triad_normalise(&self.triads[next], start);
          self.triads[tid] = triad_normalise(&self.triads[tid], start);
          log::info!("after normalise ({}): {} {} {}",
            self.triads[next].a, self.triads[next].b, self.triads[next].c,
            start,
          );
          if circumcircle_contains(&self.triads[next], self.points[outside], &self.points) {
            do_flip(tid, next, &mut self.triads, &self.points);
          }
        }
      }
    }
  }

  pub fn step(&mut self) -> std::task::Poll<Result<(), ()>> {
    //! Attempt to add the next point to the convex hull, and add triangles to any visible edges in the hull
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

    // log::info!("Hull: \n{}", self.hull.iter().map(|h| format!("{}, {}", self.points[*h].coords.x, self.points[*h].coords.y)).collect::<Vec<_>>().join("\n"));

    // visibility test: an edge (h0, h1) is visible if the angle (h0, h1, pt) is counter-clockwise.
    let visible = |hid: usize, hnext: usize, hull: &[usize], points: &[Point]| -> bool {
      let h0 = &points[hull[hid]];
      let h1 = &points[hull[hnext]];
      let angle = (h1-h0).perp(&(pt-h1));
      let visible = angle > 0f32;
      // log::info!(
      //   "Edge {}: ({}, {}) -> {}: ({}, {}) is {} ",
      //   hid, h0.coords.x, h0.coords.y,
      //   hnext, h1.coords.x, h1.coords.y,
      //   if visible { "visible" } else { "not visible" }
      // );
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

      log::info!("Found new visible edge");

      let h0 = self.hull[hid];
      let h1 = self.hull[hnext];
      
      let (centre, radius) = circumcircle(h0, h1, pid, &self.points).unwrap();
      let tid = self.triads.len();
      let mut triad = Triad::new(h0, h1, pid, centre, radius);
      self.triads.push(triad);
      log::info!("new triad with vertices {}, {}, {}", h0, h1, pid);

      for (a, b) in [(h0, h1), (h1, pid), (pid, h0)] {
        if let Some(other) = self.triads_by_edge.get(&sort_edge_lex(a, b)) {
          log::info!("new neighbour with vertices {}, {}, {}", self.triads[*other].a, self.triads[*other].b, self.triads[*other].c);
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

    let mut triads_by_edge = BTreeMap::new();
    triads_by_edge.insert(sort_edge_lex(a, b), 0);
    triads_by_edge.insert(sort_edge_lex(b, c), 0);
    triads_by_edge.insert(sort_edge_lex(c, a), 0);

    Delaunay {
      origin: points[0].clone(),
      hull: vec![ a, b, c, ],
      sorted_points: sorted,
      next: 0,
      triads_by_edge,

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
      writeln!(out, "$neighbours << EOD")?;
      for triad in &self.triads {
        for n in [triad.ab, triad.bc, triad.ca] {
          if let Some(n) = n {
            let centre = triad.com(&self.points);
            writeln!(out, "{} {} {} {}", 
              centre.coords.x, centre.coords.y,
              Point::from(self.triads[n].com(&self.points) - centre).coords.x,
              Point::from(self.triads[n].com(&self.points) - centre).coords.y
            );
          }
        }
      }
      writeln!(out, "EOD")?;
      writeln!(out, "plot $points with points, \\")?;
      writeln!(out, "$edges with vectors nohead, \\")?;
      writeln!(out, "$neighbours with vectors, \\")?;
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