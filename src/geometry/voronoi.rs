use {
  super::{
    Point,
    point_cmp_lex_yx,
    circumcircle,
  },
  std::{
    io,
    collections::{
      BTreeSet,
      VecDeque
    }
  },
};

pub struct Edge {
  src: usize, // index into vertex list
  dst: usize, // index into vertex list
}

// note: at any time we can compute the polynomial given
// its focus and the position of the scanline, which acts as a directrix
pub struct Segment {
  pub parent: Option<usize>,
  pub children: (Option<usize>, Option<usize>),
  pub focus: usize, // index into a point list
  pub left: Option<usize>,
  pub right: Option<usize>,
}

impl Segment {
  fn eqn(&self, directrix: Point, points: &[Point]) -> String {
    // parabola formula: (x-h)^{2} = 4*p*(y-k)
    // where the focus is (h, k+p) and directrix is the line y=k-p
    //
    // from the above we have:
    let focus = points[self.focus];
    let h = focus.coords.x;
    let p = (focus.coords.y - directrix.coords.y) / 2f32;
    let k = focus.coords.y - p;

    format!("((x-{})**2 / -4) + {}", h, k)
  }

  fn breakpoint_with(
    &self,
    other: &Segment,
    directrix: Point,
    points: &[Point],
  ) -> f32 {
    let f1 = &points[self.focus];
    let (x1, y1) = (f1.coords.x, f1.coords.y);
    let f2 = &points[other.focus];
    let (x2, y2) = (f2.coords.x, f2.coords.y);
    let l = directrix.coords.y;

    let d1 = 1f32 / (2f32 * (x1 - l));
    let d2 = 1f32 / (2f32 * (x2 - l));
    let a = d1 - d2;
    let b = 2f32 * (x2 * d2 - x1 * d1);
    let c = (y1 * y1 + x1 * x1 - l * l) * d1 - (y2 * y2 + x2 * x2 - l * l) * d2;
    let delta = b * b - 4f32 * a * c;
    (-b * delta.sqrt()) / (2f32 * a)
  }

  fn intercept(&self, directrix: Point, points: &[Point]) -> Point {
    // parabola formula: (x-h)^{2} = 4*p*(y-k)
    // where the focus is (h, k+p) and directrix is the line y=k-p
    //
    // from the above we have:
    let focus = points[self.focus];
    let h = focus.coords.x;
    let p = (focus.coords.y - directrix.coords.y) / 2f32;
    let k = focus.coords.y - p;

    // if we rearrange the parabola formula above to solve for y, we get
    // y = ((x-h)^{2} / 4) + k
    Point::new(
      directrix.coords.x,
      (((directrix.coords.x - h) * (directrix.coords.x - h)) / 4f32) + k,
    )
  }
}

pub struct Voronoi {
  pub points: Vec<Point>,
  pub vertices: Vec<Point>,
  pub edges: Vec<Edge>,
  pub segments: Vec<Segment>,
  pub live_segments: BTreeSet<usize>,
}

impl Voronoi {
  pub fn new(points: &[Point], max_events: usize) -> Voronoi {
    //! Returns a list of pairs:
    //!   (point, polygon)
    //! where `polygon` is the set of points defining the Voronoi polygon
    //! around `point`.

    let mut sorted = (0..points.len()).collect::<Vec<usize>>();
    sorted.sort_by(|a, b| point_cmp_lex_yx(&points[*a], &points[*b]).unwrap());

    // These are the nodes in the voronoi diagram - keep terminology
    // consistent to avoid confusion:
    //  Points -> in the original point set
    //  Vertices -> in the voronoi diagram
    let mut vertices: Vec<Point> = vec![];
    let mut edges: Vec<Edge> = vec![];
    let mut segments: Vec<Segment> = vec![];

    fn find_segment(
      directrix: Point,
      segments: &[Segment],
      points: &[Point],
    ) -> Option<usize> {
      let mut nid = 0;
      loop {
        if nid >= segments.len() {
          break None;
        }
        let segment = &segments[nid];
        let break_left = segment
          .children
          .0
          .map(|id| segments[id].breakpoint_with(&segment, directrix, points))
          .unwrap_or(f32::MIN);
        let break_right = segment
          .children
          .1
          .map(|id| segments[id].breakpoint_with(&segment, directrix, points))
          .unwrap_or(f32::MAX);
        if directrix.coords.x < break_left {
          nid = segment.children.0.unwrap();
        } else if directrix.coords.x > break_right {
          nid = segment.children.1.unwrap();
        } else {
          break Some(nid);
        }
      }
    }

    enum Event {
      Point { eid: usize, idx: usize },
      Segment { eid: usize, idx: usize },
    }

    let mut next_eid: usize = 0;
    let mut events: VecDeque<Event> = sorted
      .iter()
      .map(|i| {
        next_eid += 1;
        Event::Point {
          eid: next_eid - 1,
          idx: *i,
        }
      })
      .collect::<VecDeque<_>>();
    let mut live_segments: BTreeSet<usize> = BTreeSet::new();

    while let Some(event) = events.pop_front() {
      match event {
        Event::Point { eid, idx } => {
          if eid > max_events {
            break;
          }
          // Find the segment directly above point `idx`
          // We could speed this up a lot with a smarter data structure
          // but implementing that data structure is out of scope for now
          let point = points[idx];
          log::info!("Point event for ({}, {})", point.coords.x, point.coords.y);
          match find_segment(point, &segments, &points) {
            None => {
              log::info!("First segment");
              segments.push(Segment {
                parent: None,
                children: (None, None),
                focus: idx,
                left: None,
                right: None,
              });
              live_segments.insert(segments.len() - 1);
            }
            Some(segid) => {
              // this segment will no longer be live, we're about to split it
              live_segments.remove(&segid);

              // two new vertices: currently these are at the same point but they'll separate
              // later on
              let a0p = segments[segid].intercept(points[idx], &points);
              let a1p = a0p;
              let a0 = vertices.len();
              vertices.push(a0p);
              let a1 = vertices.len();
              vertices.push(a1p);

              edges.push(Edge { src: a0, dst: a1 });
              edges.push(Edge { src: a1, dst: a0 });

              let lhs_id = segments.len();
              let mid_id = lhs_id + 1;
              let rhs_id = lhs_id + 2;

              let lhs = Segment {
                parent: Some(segid),
                children: (None, None),
                focus: segments[segid].focus,
                left: segments[segid].left,
                right: Some(mid_id),
              };
              let mid = Segment {
                parent: Some(segid),
                children: (None, None),
                focus: idx,
                left: Some(lhs_id),
                right: Some(rhs_id),
              };
              let rhs = Segment {
                parent: Some(segid),
                children: (None, None),
                focus: segments[segid].focus,
                left: Some(mid_id),
                right: segments[segid].right,
              };
              for seg in [lhs, mid, rhs] {
                segments.push(seg);
                live_segments.insert(segments.len() - 1);
              }
              if let Some(llhs) = segments[segid].left {
                log::info!("Found llhs, generating segment event");
                segments[llhs].right = Some(lhs_id);
                // Check foci of llhs, lhs, mid for circle event
                let (a, b, c) = (
                  segments[llhs].focus,
                  segments[lhs_id].focus,
                  segments[mid_id].focus,
                );
                let (_, y) = circumcircle(a, b, c, &points).unwrap();
                events.push_front(Event::Segment {
                  eid: next_eid,
                  idx: lhs_id,
                });
                next_eid += 1;
              }
              if let Some(rrhs) = segments[segid].right {
                log::info!("Found rrhs, generating segment event");
                segments[rrhs].left = Some(rhs_id);
                // check foci of rrhs, rhs, mid for circle event
                let (a, b, c) = (
                  segments[mid_id].focus,
                  segments[rhs_id].focus,
                  segments[rrhs].focus,
                );
                let (_, y) = circumcircle(a, b, c, &points).unwrap();
                events.push_front(Event::Segment {
                  eid: next_eid,
                  idx: rhs_id,
                });
                next_eid += 1;
              }
            }
          }
        }
        Event::Segment { eid, idx } => {
          log::info!("Segment event for {}", idx);
          live_segments.remove(&idx);
          // update breakpoints?
          for segid in &live_segments {}
        }
      }
    }

    Voronoi {
      points: points.to_vec(),
      vertices,
      edges,
      segments,
      live_segments,
    }
  }

  pub fn gnuplot(&self, directrix: Point) -> io::Result<String> {
    use io::{
      Read,
      Seek,
      Write,
    };
    let mut out = io::Cursor::new(Vec::new());
    writeln!(out, "set xrange [0:10]")?;
    writeln!(out, "set yrange [0:10]")?;
    writeln!(out, "$points << EOD")?;
    for point in &self.points {
      writeln!(out, "{} {} 2", point.coords.x, point.coords.y)?;
    }
    writeln!(out, "EOD")?;

    writeln!(out, "$vertices << EOD")?;
    for edge in &self.edges {
      writeln!(
        out,
        "{} {} {} {}",
        self.vertices[edge.src].coords.x,
        self.vertices[edge.src].coords.y,
        self.vertices[edge.dst].coords.x - self.vertices[edge.src].coords.x,
        self.vertices[edge.dst].coords.y - self.vertices[edge.src].coords.y
      )?;
    }
    writeln!(out, "EOD")?;

    writeln!(out, "set yrange [*:*] reverse")?;
    writeln!(out, "plot $points pointtype 7 with points, \\")?;

    for segment in &self.live_segments {
      let segment = &self.segments[*segment];
      writeln!(
        out,
        "{} title \"focus({}, {})\" with lines, \\",
        segment.eqn(directrix, &self.points),
        self.points[segment.focus].coords.x,
        self.points[segment.focus].coords.y
      )?;
    }

    writeln!(out, "$vertices with vectors,\\")?;

    writeln!(
      out,
      "{} title \"directrix({}, {})\"\\",
      directrix.coords.y, directrix.coords.x, directrix.coords.y
    )?;

    writeln!(out);

    let mut result = String::new();
    out.seek(io::SeekFrom::Start(0))?;
    out.read_to_string(&mut result)?;
    Ok(result)
  }
}

#[cfg(test)]
#[test]
fn parabolas() {
  let points = &[
    Point::new(2f32, 3f32), // focus
  ];
  let segment = Segment {
    // note: breakpoints aren't used in this test, and these values are invalid.
    focus: 0,
    left: None,
    right: None,
    children: (None, None),
    parent: None,
  };

  assert_eq!(
    segment.intercept(Point::new(2f32, 9f32), points),
    Point::new(2f32, 6f32)
  );
  assert_eq!(
    segment.intercept(Point::new(4f32, 9f32), points),
    Point::new(4f32, 7f32)
  );
}