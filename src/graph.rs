use crate::idalloc::IdAlloc;
use std::{
  collections::{
    btree_map::Entry,
    BTreeMap,
    BTreeSet,
  },
  io::{
    self,
    Write,
  },
};

pub struct Node<T> {
  ingress: BTreeSet<usize>,
  egress: BTreeSet<usize>,
  value: T,
}

impl<T> Node<T> {
  fn new(value: T) -> Node<T> {
    Node {
      ingress: Default::default(),
      egress: Default::default(),
      value,
    }
  }
}

pub struct Edge<T> {
  src: usize,
  dest: usize,
  value: T,
}

/// This directed graph is constructed as two mappings:
///   Node ID -> Node { IDs of incoming edges,
///                     IDs of outgoing edges }
///   Edge ID -> Edge { ID of source node,
///                     ID of destination node }
///
/// There are a couple of obvious extension points:
///   - Adding data to the nodes
///   - Adding data (e.g. costs / weights) to edges
#[derive(Default)]
pub struct Graph<N, E> {
  nodes: BTreeMap<usize, Node<N>>,
  edges: BTreeMap<usize, Edge<E>>,
  id_alloc: IdAlloc,
}

impl<N> Graph<N, ()> {
  pub fn add_edges(&mut self, src: usize, dests: &[usize]) -> Vec<usize> {
    dests
      .iter()
      .filter_map(|dest| self.add_edge(src, *dest, ()))
      .collect()
  }
}

impl<N, E> Graph<N, E> {
  pub fn get_node(&self, nid: usize) -> Option<&N> {
    self.nodes.get(&nid).map(|n| &n.value)
  }

  // Returns Some(nid) on success, None if the node already exists
  // (ideally we'd allocate `nid` ourselves and reuse them, but I'm trying to keep
  // this simple)
  pub fn add_node(&mut self, value: N) -> Option<usize> {
    let nid = self.id_alloc.alloc_id();
    match self.nodes.entry(nid) {
      Entry::Vacant(entry) => {
        entry.insert(Node::new(value));
        Some(nid)
      }
      Entry::Occupied(_) => {
        panic!("Allocated an already-in-use id?");
      }
    }
  }

  pub fn remove_node(&mut self, nid: usize) -> Option<usize> {
    match self.nodes.entry(nid) {
      Entry::Vacant(_) => panic!("double-free of id {}", nid),
      Entry::Occupied(entry) => {
        self.id_alloc.free_id(nid);
        let node = entry.remove();
        for edge in node.ingress.into_iter().chain(node.egress.into_iter()) {
          self.edges.remove(&edge);
        }
        Some(nid)
      }
    }
  }

  pub fn add_edge(&mut self, src: usize, dest: usize, value: E) -> Option<usize> {
    if !self.nodes.contains_key(&src) || !self.nodes.contains_key(&dest) {
      return None;
    }

    let eid = self.id_alloc.alloc_id();
    match self.edges.entry(eid) {
      Entry::Occupied(_) => {
        panic!("Allocated an already-in-use id?");
      }
      Entry::Vacant(entry) => {
        entry.insert(Edge { src, dest, value });
        self.nodes.get_mut(&src).unwrap().egress.insert(eid);
        self.nodes.get_mut(&dest).unwrap().ingress.insert(eid);
        return Some(eid);
      }
    }
  }

  pub fn remove_edge(&mut self, eid: usize) -> Option<usize> {
    match self.edges.entry(eid) {
      Entry::Vacant(_) => panic!("double-free of id {}", eid),
      Entry::Occupied(entry) => {
        self
          .nodes
          .get_mut(&entry.get().src)
          .unwrap()
          .egress
          .remove(&eid);
        self
          .nodes
          .get_mut(&entry.get().dest)
          .unwrap()
          .ingress
          .remove(&eid);
        entry.remove();
        self.id_alloc.free_id(eid);
        Some(eid)
      }
    }
  }

  pub fn djikstra_shortest_path(&self, source: usize, dest: usize) -> Option<Vec<usize>> {
    use std::cmp::Reverse;

    let mut dists: BTreeMap<usize, (usize, usize)> = BTreeMap::new(); // Node ID -> Distance, Predecessor
    dists.entry(source).or_insert((0, source));

    let mut q = self.nodes.keys().cloned().collect::<Vec<_>>();
    q.sort_by_key(|nid| Reverse(dists.get(nid).map(|t| t.0).unwrap_or(usize::MAX)));

    while let Some(nid) = q.pop() {
      let cost = match dists.get(&nid) {
        None => break,
        Some((cost, _)) => *cost,
      };
      for edge in self
        .nodes
        .get(&nid)
        .unwrap()
        .egress
        .iter()
        .map(|eid| self.edges.get(eid).unwrap())
      {
        let cost = cost + 1;
        if cost < dists.get(&edge.dest).map(|t| t.0).unwrap_or(usize::MAX) {
          dists.insert(edge.dest, (cost, nid));
        }
      }

      q.sort_by_key(|nid| Reverse(dists.get(&nid).map(|t| t.0).unwrap_or(usize::MAX)));
    }

    let mut result = vec![dest];
    let mut n = dest;
    while n != source {
      n = dists.get(&n).unwrap().1;
      result.push(n);
    }

    Some(result.into_iter().rev().collect())
  }

  pub fn is_tree(&self) -> bool {
    // Algorithm:
    // - Start from a random node
    // - Verify only one ingress node (if 0, identify the root of the tree)
    // - FOR ingress edge:
    //    - Add source to the explore list, if it hasn't been explored
    // - FOR egress edges:
    //    - Remove destination from node set. If it's not there, it's already been
    //      visited by an egress edge and the graph is not  a tree. (or it's the
    //      start node)
    //    - Add the destination to the explore list.
    // - Pick a new edge from the explore list.

    let mut nodes = self.nodes.keys().collect::<BTreeSet<_>>();
    let mut explore = Vec::new();

    // An empty graph is a tree by definition
    let start = if let Some(nid) = nodes.pop_first() {
      nid
    } else {
      return true;
    };
    explore.push(start);

    while let Some(nid) = explore.pop() {
      let node = &self.nodes.get(nid).unwrap();
      match node.ingress.len() {
        0 => {
          // we found the root of the tree
        }
        1 => {
          let edge = self.edges.get(node.ingress.first().unwrap()).unwrap();
          if nodes.contains(&edge.src) {
            nodes.remove(&edge.src);
            explore.push(&edge.src);
          }
        }
        _ => return false,
      }

      for eid in node.egress.iter() {
        let edge = &self.edges[eid];
        match nodes.contains(&edge.dest) {
          false if edge.dest == *start => continue,
          false => return false,
          true => {
            nodes.remove(&edge.dest);
            explore.push(&edge.dest);
          }
        }
      }
    }

    return nodes.is_empty();
  }

  pub fn is_eulerian(&self) -> bool {
    //! Requirements
    //! - Graph is connected
    //! - All nodes have even indegree
    //!
    //! If the graph didn't have to be connected, we could just check the
    //! indegree of every node directly Instead we have to walk the graph to
    //! make sure it's connected.
    //!
    //! Algorithm:
    //! - Walk the graph starting from a random node
    //! - When a node is visited, remove it from the set of unvisited nodes
    //!   Visit each node at most once
    //! - If any nodes are unvisited at the end, return false
    //! - If any node has an odd ingress count, return false.

    let mut nodes = self.nodes.keys().collect::<BTreeSet<_>>();

    let mut explore = Vec::new();
    let start = nodes.pop_first().unwrap();
    explore.push(*start);

    while let Some(nid) = explore.pop() {
      let node = &self.nodes.get(&nid).unwrap();
      for eid in node.egress.iter() {
        let edge = &self.edges.get(&eid).unwrap();
        if nodes.remove(&edge.dest) {
          explore.push(edge.dest);
        }
      }
    }

    return nodes.is_empty()
      && self
        .nodes
        .keys()
        .all(|nid| self.nodes.get(nid).unwrap().ingress.len() % 2 == 0);
  }

  pub fn write_adjacency_matrix(&self, writer: &mut dyn Write) -> io::Result<()> {
    write!(writer, "|   |")?;
    for i in self.nodes.keys() {
      write!(writer, " {} |", i)?;
    }
    writeln!(writer)?;
    for i in self.nodes.keys() {
      write!(writer, "| {} |", i)?;
      for j in self.nodes.keys() {
        if self
          .nodes
          .get(i)
          .unwrap()
          .egress
          .iter()
          .map(|eid| self.edges.get(eid).unwrap().dest)
          .any(|it| it == *j)
        {
          write!(writer, " x |")?;
        } else {
          write!(writer, " _ |")?;
        }
      }
      writeln!(writer)?;
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::Graph;
  use std::io;
  #[test]
  fn identifies_trees() {
    let mut g = Graph::default();

    let nids = (0..7)
      .into_iter()
      .filter_map(|i| g.add_node(()))
      .collect::<Vec<_>>();

    for i in (0usize..7) {
      assert!(nids.contains(&i));
    }

    let eids = &[(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
      .into_iter()
      .map(|(s, d)| g.add_edge(*s, *d, ()))
      .collect::<Vec<_>>();

    assert!(g.is_tree());

    // add a cycle to the tree
    let cycle = g.add_edge(4, 0, ()).unwrap();
    assert_eq!(g.is_tree(), false);
    g.remove_edge(cycle);
    assert_eq!(g.is_tree(), true);

    // add an extra node to the tree
    g.add_node(());
    assert_eq!(g.is_tree(), false);
  }

  #[test]
  fn identifies_eulerian_circuits() {
    let mut g = Graph::default();

    for i in 0..7 {
      g.add_node(());
    }

    g.add_edges(0, &[1, 2]);
    g.add_edges(1, &[0, 2, 3, 4]);
    g.add_edges(2, &[0, 1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5, 6]);
    g.add_edges(5, &[2, 3, 4, 6]);
    g.add_edges(6, &[4, 5]);

    assert!(g.is_eulerian());
    assert_eq!(g.is_tree(), false);

    // ---

    let mut g = Graph::default();
    for i in 0..6 {
      g.add_node(());
    }

    g.add_edges(0, &[1, 2]);
    g.add_edges(1, &[0, 2, 3, 4]);
    g.add_edges(2, &[0, 1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5]);
    g.add_edges(5, &[2, 3, 4]);

    assert_eq!(g.is_eulerian(), false);
    assert_eq!(g.is_tree(), false);

    // ---

    let mut g = Graph::default();
    for i in 0..6 {
      g.add_node(());
    }

    g.add_edges(1, &[2, 3, 4]);
    g.add_edges(2, &[1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5]);
    g.add_edges(5, &[2, 3, 4]);

    assert_eq!(g.is_eulerian(), false);
    assert_eq!(g.is_tree(), false);

    // --- this is Euler's Konigsberg problem

    for i in 0..4 {
      g.add_node(());
    }
    g.add_edges(0, &[1, 3]);
    g.add_edges(1, &[0, 2, 3]);
    g.add_edges(2, &[1, 3]);
    g.add_edges(3, &[0, 1, 2]);

    assert_eq!(g.is_eulerian(), false);
    assert_eq!(g.is_tree(), false);
  }

  fn verify_matrix_printout<N, E>(graph: &mut Graph<N, E>, expectation: &str) -> io::Result<()> {
    let mut buffer = Vec::new();
    graph.write_adjacency_matrix(&mut buffer)?;
    let matrix = std::io::read_to_string(&mut buffer.as_slice())?;
    assert_eq!(matrix, expectation);

    Ok(())
  }

  #[test]
  fn prints_adjacency_matrices() -> io::Result<()> {
    let mut g = Graph::default();

    for i in 0..7 {
      g.add_node(());
    }

    g.add_edges(0, &[1, 2]);
    g.add_edges(1, &[0, 2, 3, 4]);
    g.add_edges(2, &[0, 1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5, 6]);
    g.add_edges(5, &[2, 3, 4, 6]);
    g.add_edges(6, &[4, 5]);

    verify_matrix_printout(
      &mut g,
      "\
      |   | 0 | 1 | 2 | 3 | 4 | 5 | 6 |\n| 0 | _ | x | x | _ | _ | _ | _ |\n| 1 | x | _ | x | x | \
       x | _ | _ |\n| 2 | x | x | _ | x | _ | x | _ |\n| 3 | _ | x | x | _ | x | x | _ |\n| 4 | _ \
       | x | _ | x | _ | x | x |\n| 5 | _ | _ | x | x | x | _ | x |\n| 6 | _ | _ | _ | _ | x | x \
       | _ |\n",
    )?;

    // ---

    let mut g: Graph<(), ()> = Graph::default();
    for i in 0..6 {
      g.add_node(());
    }

    g.add_edges(0, &[1, 2]);
    g.add_edges(1, &[0, 2, 3, 4]);
    g.add_edges(2, &[0, 1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5]);
    g.add_edges(5, &[2, 3, 4]);

    verify_matrix_printout(
      &mut g,
      "\
      |   | 0 | 1 | 2 | 3 | 4 | 5 |\n| 0 | _ | x | x | _ | _ | _ |\n| 1 | x | _ | x | x | x | _ \
       |\n| 2 | x | x | _ | x | _ | x |\n| 3 | _ | x | x | _ | x | x |\n| 4 | _ | x | _ | x | _ | \
       x |\n| 5 | _ | _ | x | x | x | _ |\n",
    )?;

    // ---

    let mut g = Graph::default();
    for i in 0..6 {
      g.add_node(());
    }

    g.add_edges(1, &[2, 3, 4]);
    g.add_edges(2, &[1, 3, 5]);
    g.add_edges(3, &[1, 2, 4, 5]);
    g.add_edges(4, &[1, 3, 5]);
    g.add_edges(5, &[2, 3, 4]);

    verify_matrix_printout(
      &mut g,
      "\
      |   | 0 | 1 | 2 | 3 | 4 | 5 |\n| 0 | _ | _ | _ | _ | _ | _ |\n| 1 | _ | _ | x | x | x | _ \
       |\n| 2 | _ | x | _ | x | _ | x |\n| 3 | _ | x | x | _ | x | x |\n| 4 | _ | x | _ | x | _ | \
       x |\n| 5 | _ | _ | x | x | x | _ |\n",
    )?;

    // --- this is Euler's Konigsberg problem

    let mut g = Graph::default();

    for i in 0..4 {
      g.add_node(());
    }
    g.add_edges(0, &[1, 3]);
    g.add_edges(1, &[0, 2, 3]);
    g.add_edges(2, &[1, 3]);
    g.add_edges(3, &[0, 1, 2]);

    assert_eq!(g.is_eulerian(), false);
    assert_eq!(g.is_tree(), false);

    verify_matrix_printout(
      &mut g,
      "\
      |   | 0 | 1 | 2 | 3 |\n| 0 | _ | x | _ | x |\n| 1 | x | _ | x | x |\n| 2 | _ | x | _ | x \
       |\n| 3 | x | x | x | _ |\n",
    )?;

    Ok(())
  }

  #[test]
  fn data_example() -> io::Result<()> {
    let city_names = &[
      // 0
      "Canberra",
      // 1
      "Melbourne",
      // 2
      "Hobart",
      // 3
      "Adelaide",
      // 4
      "Perth",
      // 5
      "Darwin",
      // 6
      "Brisbane",
      // 7
      "Sydney",
      // 8
      "Auckland",
      // 9
      "Wellington",
      // 10
      "Christchurch",
      // 11
      "Dunedin",
    ];

    let mut g = Graph::default();
    for i in 0..(city_names.len()) {
      g.add_node(());
    }

    g.add_edges(0, &[1, 2, 3, 4, 5, 6, 7, 9]);
    g.add_edges(1, &[7, 9]);
    g.add_edges(7, &[1, 8]);
    g.add_edges(9, &[8, 10, 11]);
    g.add_edges(10, &[11]);

    g.write_adjacency_matrix(&mut io::stdout().lock()).unwrap();

    verify_matrix_printout(
      &mut g,
      "|   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |\n| 0 | _ | x | x | x | x | x | x | \
       x | _ | x | _ | _ |\n| 1 | _ | _ | _ | _ | _ | _ | _ | x | _ | x | _ | _ |\n| 2 | _ | _ | \
       _ | _ | _ | _ | _ | _ | _ | _ | _ | _ |\n| 3 | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | \
       _ |\n| 4 | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ |\n| 5 | _ | _ | _ | _ | _ | _ | \
       _ | _ | _ | _ | _ | _ |\n| 6 | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ |\n| 7 | _ | \
       x | _ | _ | _ | _ | _ | _ | x | _ | _ | _ |\n| 8 | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | \
       _ | _ |\n| 9 | _ | _ | _ | _ | _ | _ | _ | _ | x | _ | x | x |\n| 10 | _ | _ | _ | _ | _ | \
       _ | _ | _ | _ | _ | _ | x |\n| 11 | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ | _ |\n",
    )
    .unwrap();

    assert_eq!(g.djikstra_shortest_path(1, 11), Some(vec![1, 9, 11]));

    Ok(())
  }
}
