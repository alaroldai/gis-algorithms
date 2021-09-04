use gis_algorithms::geometry::*;
use nalgebra::Point2;
use rand::{
  prelude::SliceRandom,
  thread_rng,
  Rng,
};
use std::{
  fs,
  io::{
    self,
    Write,
  },
  path::PathBuf,
};
use structopt::StructOpt;

#[derive(StructOpt)]
#[structopt(name = "GIS Algorithms")]
struct App {
  #[structopt(short, long, parse(from_os_str), default_value = "gen")]
  output: PathBuf,

  #[structopt(long)]
  examples: usize,

  #[structopt(long)]
  delaunay_steps: bool,
}

impl App {
  fn gen_bounding_boxes(&self) -> io::Result<()> {
    let mut points = Vec::new();
    for x in 0..=10 {
      for y in 0..=10 {
        points.push(Point2::new(x as f32, y as f32));
      }
    }

    for idx in 0..self.examples {
      let sample = points
        .choose_multiple(&mut thread_rng(), points.len() / 10)
        .cloned()
        .collect::<Vec<_>>();
      let mut of = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(self.output.join(&format!("hulls-{}.dot", idx)))?;
      writeln!(of, "graph {{")?;
      writeln!(of, "  graph [layout=fdp]")?;
      for (idx, point) in sample.iter().enumerate() {
        writeln!(
          of,
          "  {} [shape=point,pos=\"{},{}!\"]",
          idx, point.x, point.y
        )?;
      }

      let hull = convex_hull(&sample);
      for i in 1..hull.len() {
        writeln!(of, "  {} -- {}", hull[i - 1], hull[i])?;
      }
      writeln!(of, "  {} -- {}", hull[hull.len() - 1], hull[0])?;
      writeln!(of, "}}")?;
    }

    Ok(())
  }

  fn gen_delaunay(&self) -> io::Result<()> {
    let mut points = Vec::new();
    for x in 0..10 {
      for y in 0..10 {
        points.push(Point2::new(
          thread_rng().gen::<f32>() + x as f32,
          thread_rng().gen::<f32>() + y as f32,
        ));
      }
    }

    for i in 0..self.examples {
      let sample = points
        .choose_multiple(&mut thread_rng(), points.len() / 10)
        .cloned()
        .collect::<Vec<_>>();

      let mut delaunay = Delaunay::new(&sample);

      for step in 0.. {
        if self.delaunay_steps {
          let mut of = fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(self.output.join(&format!("delaunay-{}-{}.plt", i, step)))?;
          write!(of, "{}", delaunay.gnuplot()?)?;
        }

        match delaunay.step() {
          std::task::Poll::Ready(Ok(())) => break,
          std::task::Poll::Ready(Err(())) => break,
          std::task::Poll::Pending => continue,
        };
      }

      let mut of = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(self.output.join(&format!("delaunay-{}.plt", i)))?;
      write!(of, "{}", delaunay.gnuplot()?)?;
    }

    Ok(())
  }

  #[cfg(feature = "voronoi")]
  fn gen_voronoi(&self) -> io::Result<()> {
    let mut points = Vec::new();
    for x in 0..10 {
      for y in 0..10 {
        points.push(Point2::new(
          thread_rng().gen::<f32>() + x as f32,
          thread_rng().gen::<f32>() + y as f32,
        ));
      }
    }

    for i in 0..self.examples {
      let mut sample = points
        .choose_multiple(&mut thread_rng(), points.len() / 10)
        .cloned()
        .collect::<Vec<_>>();

      sample.sort_by(|a, b| point_cmp_lex_yx(a, b).unwrap());

      // let voronoi = Voronoi {
      //   points: sample.clone(),
      //   vertices: vec![],
      //   edges: vec![],
      //   segments: vec![
      //     Segment {
      //       break_left: None,
      //       break_right: None,
      //       focus: 0,
      //     }
      //   ]
      // };

      let voronoi = Voronoi::new(&sample, 3);

      let mut of = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(self.output.join(&format!("voronoi-{}.plt", i)))?;
      write!(of, "{}", voronoi.gnuplot(sample[3])?)?;
    }
    Ok(())
  }
}

fn main() -> io::Result<()> {
  env_logger::init();

  let app = App::from_args();

  fs::create_dir_all(&app.output)?;

  app.gen_bounding_boxes()?;

  #[cfg(feature = "voronoi")]
  {
    app.gen_voronoi()?;
  }

  app.gen_delaunay()?;

  Ok(())
}
