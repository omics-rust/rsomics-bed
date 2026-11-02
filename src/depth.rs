use std::collections::BTreeMap;

#[derive(Clone, Copy)]
pub(crate) struct DepthRun {
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) depth: usize,
}

#[derive(Default)]
pub(crate) struct DepthBuffer {
    events: Vec<(u64, bool)>,
    runs: Vec<DepthRun>,
    histogram: BTreeMap<usize, u64>,
}

impl DepthBuffer {
    pub(crate) fn clear_events(&mut self) {
        self.events.clear();
    }

    pub(crate) fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.events.push((start, true));
            self.events.push((end, false));
        }
    }

    pub(crate) fn reduce(&mut self, start: u64, end: u64) {
        self.runs.clear();
        self.histogram.clear();
        self.events.sort_unstable_by_key(|event| event.0);
        let mut cursor = start;
        let mut depth = 0;
        let mut event = 0;
        while event < self.events.len() {
            let position = self.events[event].0;
            self.push_run(cursor, position, depth);
            let mut starts = 0;
            let mut ends = 0;
            while event < self.events.len() && self.events[event].0 == position {
                if self.events[event].1 {
                    starts += 1;
                } else {
                    ends += 1;
                }
                event += 1;
            }
            depth = depth + starts - ends;
            cursor = position;
        }
        self.push_run(cursor, end, depth);
    }

    pub(crate) fn covered_bases(&self, length: u64) -> u64 {
        length - self.histogram.get(&0).copied().unwrap_or(0)
    }

    pub(crate) fn runs(&self) -> &[DepthRun] {
        &self.runs
    }

    pub(crate) fn histogram(&self) -> &BTreeMap<usize, u64> {
        &self.histogram
    }

    fn push_run(&mut self, start: u64, end: u64, depth: usize) {
        if start < end {
            self.runs.push(DepthRun { start, end, depth });
            *self.histogram.entry(depth).or_default() += end - start;
        }
    }
}
