use log::*;
use std::{
  cmp::Ordering,
  collections::VecDeque,
  fmt::{
    self,
    Debug,
  },
};

#[derive(PartialOrd, PartialEq)]
struct RangeInclusive {
  pub start: usize,
  pub end: usize,
}

impl Debug for RangeInclusive {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
    write!(f, "[{}, {}]", self.start, self.end)
  }
}

impl Iterator for RangeInclusive {
  type Item = usize;
  fn next(&mut self) -> Option<usize> {
    if self.start <= self.end {
      self.start += 1;
      return Some(self.start - 1);
    }
    return None;
  }
}

trait SliceExt<T> {
  fn replace(&mut self, idx: usize, expand: impl FnOnce(Option<T>) -> Vec<T>);
}

impl<T> SliceExt<T> for VecDeque<T> {
  fn replace(&mut self, idx: usize, expand: impl FnOnce(Option<T>) -> Vec<T>) {
    let ins = expand(self.remove(idx));
    for it in ins.into_iter().rev() {
      self.insert(idx, it);
    }
  }
}

pub struct IdAlloc {
  // list of free ranges, less than alloc_limit
  alloc_free_range: VecDeque<RangeInclusive>,
}

impl Default for IdAlloc {
  fn default() -> IdAlloc {
    IdAlloc {
      alloc_free_range: VecDeque::from(vec![RangeInclusive {
        start: 0,
        end: usize::MAX,
      }]),
    }
  }
}

impl IdAlloc {
  pub fn alloc_id(&mut self) -> usize {
    match self.alloc_free_range.pop_front() {
      None => {
        panic!("out of identifiers");
      }
      Some(mut front_range) => {
        let result = front_range.next().unwrap();
        self.alloc_free_range.push_front(front_range);

        log::trace!("alloc {}", result);

        return result;
      }
    }
  }

  pub fn ensure_id(&mut self, id: usize) {
    match self.alloc_free_range.binary_search_by(|a| {
      if a.start > id {
        return Ordering::Greater;
      }
      if a.end <= id {
        return Ordering::Less;
      }
      return Ordering::Equal;
    }) {
      Ok(ipt) => {
        // id is in a free list, we need to split this node
        self.alloc_free_range.replace(ipt, |existing| {
          match existing {
            Some(existing) => {
              vec![
                RangeInclusive {
                  start: existing.start,
                  end: id - 1,
                },
                RangeInclusive {
                  start: id + 1,
                  end: existing.end,
                },
              ]
            }
            None => panic!("??"),
          }
        });
      }
      Err(ipt) => {
        // id is already allocated, nothing to do
      }
    }
  }

  pub fn free_id(&mut self, id: usize) {
    // find a free range that should contain the id, if one exists

    match self.alloc_free_range.binary_search_by(|a| {
      if a.start > id {
        return Ordering::Greater;
      }
      if a.end <= id {
        return Ordering::Less;
      }
      return Ordering::Equal;
    }) {
      Ok(_) => {
        panic!("id has already been freed");
      }
      Err(ipt) => {
        let mut ins = RangeInclusive { start: id, end: id };
        // self.alloc_free_range.insert(ipt, id..=id);
        if ipt < self.alloc_free_range.len() && self.alloc_free_range[ipt].start - 1 == id {
          ins.end = self.alloc_free_range[ipt].end;
          self.alloc_free_range.remove(ipt);
        }
        if ipt > 0 && self.alloc_free_range[ipt - 1].end == id - 1 {
          ins.start = self.alloc_free_range[ipt - 1].start;
          self.alloc_free_range.insert(ipt, ins);
          self.alloc_free_range.remove(ipt - 1);
        } else {
          self.alloc_free_range.insert(ipt, ins);
        }
      }
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  #[test]
  fn it_works() {
    let mut alloc = IdAlloc::default();
    for i in 0..10 {
      assert_eq!(alloc.alloc_id(), i);
    }

    for i in 0..10 {
      alloc.free_id(i);
    }
    assert_eq!(
      alloc.alloc_free_range,
      &[RangeInclusive {
        start: 0,
        end: usize::MAX
      }]
    );

    for i in 0..10 {
      assert_eq!(alloc.alloc_id(), i);
    }

    for i in &[2usize, 4, 6, 8, 1, 3, 5, 7, 9, 0] {
      alloc.free_id(*i);
    }
    assert_eq!(
      alloc.alloc_free_range,
      &[RangeInclusive {
        start: 0,
        end: usize::MAX
      }]
    );
  }
}
