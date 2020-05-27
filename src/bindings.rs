#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

pub use super::generated_bindings::*;
use super::generated_bindings;
use num_integer::Integer;

// bindgen doesn't get the types right for the following constants
pub const FALSE: boolean = generated_bindings::FALSE as boolean;
pub const TRUE: boolean = generated_bindings::TRUE as boolean;
pub const CONSOLWIDTH: i32 = generated_bindings::CONSOLWIDTH as i32;

// bindgen gets this wrong somehow? linker won't find it.
#[allow(non_upper_case_globals)]
pub const bit: [set; 64] = [
    0o1000000000000000000000,0o400000000000000000000,
    0o200000000000000000000,0o100000000000000000000,
    0o40000000000000000000,0o20000000000000000000,
    0o10000000000000000000,0o4000000000000000000,
    0o2000000000000000000,0o1000000000000000000,
    0o400000000000000000,0o200000000000000000,
    0o100000000000000000,0o40000000000000000,0o20000000000000000,
    0o10000000000000000,0o4000000000000000,0o2000000000000000,
    0o1000000000000000,0o400000000000000,0o200000000000000,
    0o100000000000000,0o40000000000000,0o20000000000000,
    0o10000000000000,0o4000000000000,0o2000000000000,
    0o1000000000000,0o400000000000,0o200000000000,0o100000000000,
    0o40000000000,0o20000000000,0o10000000000,0o4000000000,
    0o2000000000,0o1000000000,0o400000000,0o200000000,0o100000000,
    0o40000000,0o20000000,0o10000000,0o4000000,0o2000000,0o1000000,
    0o400000,0o200000,0o100000,0o40000,0o20000,0o10000,0o4000,
    0o2000,0o1000,0o400,0o200,0o100,0o40,0o20,0o10,0o4,0o2,0o1
];

pub fn SETWORDSNEEDED(n: usize) -> usize {
    n.div_ceil(&(WORDSIZE as usize))
}

impl std::default::Default for optionblk {
    fn default() -> Self {
        optionblk{
            getcanon: 0,
            digraph: FALSE,
            writeautoms: FALSE,
            writemarkers: FALSE,
            defaultptn: TRUE,
            cartesian: FALSE,
            linelength: CONSOLWIDTH,
            outfile: std::ptr::null_mut(),
            userrefproc: None,
            userautomproc: None,
            userlevelproc: None,
            usernodeproc: None,
            usercanonproc: None,
            invarproc: None,
            tc_level: 100,
            mininvarlevel: 0,
            maxinvarlevel: 1,
            invararg: 0,
            dispatch: unsafe {&mut dispatch_graph},
            schreier: FALSE,
            extra_options: std::ptr::null_mut(),
        }
    }
}

impl std::default::Default for TracesOptions {
    fn default() -> Self {
        TracesOptions{
            getcanon: FALSE,
            writeautoms: FALSE,
            cartesian: FALSE,
            digraph: FALSE,
            defaultptn: TRUE,
            linelength: 0,
            outfile: std::ptr::null_mut(),
            strategy: 0,
            verbosity: 0,
            generators: std::ptr::null_mut(),
            userautomproc: None,
            usercanonproc: None,
            weighted: FALSE,
        }
    }
}

impl std::default::Default for statsblk {
    fn default() -> Self {
        statsblk{
            canupdates:     Default::default(),
            errstatus:      Default::default(),
            grpsize1:       Default::default(),
            grpsize2:       Default::default(),
            invapplics:     Default::default(),
            invarsuclevel:  Default::default(),
            invsuccesses:   Default::default(),
            maxlevel:       Default::default(),
            numbadleaves:   Default::default(),
            numgenerators:  Default::default(),
            numnodes:       Default::default(),
            numorbits:      Default::default(),
            tctotal:        Default::default(),
        }
    }
}

impl std::default::Default for TracesStats {
    fn default() -> Self {
        TracesStats{
            grpsize1:       Default::default(),
            grpsize2:       Default::default(),
            numgenerators:  Default::default(),
            numorbits:      Default::default(),
            treedepth:      Default::default(),
            canupdates:     Default::default(),
            errstatus:      Default::default(),
            numnodes:       Default::default(),
            interrupted:    Default::default(),
            peaknodes:      Default::default(),
        }
    }
}

pub fn empty_graph(m: usize, n: usize) -> Vec<graph> {
    vec![0; m*n]
}

pub fn uninit_graph(m: usize, n: usize) -> Vec<graph> {
    Vec::with_capacity(m*n)
}

pub fn ADDONEEDGE(g: &mut [graph], v: usize, w: usize, m: usize) {
    ADDONEARC(g, v, w, m);
    ADDONEARC(g, w, v, m);
}

pub fn ADDONEARC(g: &mut [graph], v: usize, w: usize, m: usize) {
    ADDELEMENT(GRAPHROW(g as &mut [set], v, m), w);
}

pub fn GRAPHROW(g: &mut [set], v: usize, m: usize) -> &mut [set] {
    &mut g[m * v..]
}

pub fn ADDELEMENT(setadd: &mut [set], pos: usize) {
    setadd[SETWD(pos)] |= bit[SETBT(pos)]
}

pub fn SETWD(pos: usize) -> usize {
    pos / (WORDSIZE as usize)
}

pub fn SETBT(pos: usize) -> usize {
    pos & (WORDSIZE as usize - 1)
}

impl std::default::Default for sparsegraph {
    fn default() -> Self {
        sparsegraph{
            nde:  Default::default(),
            v:    std::ptr::null_mut(),
            nv:   Default::default(),
            d:    std::ptr::null_mut(),
            e:    std::ptr::null_mut(),
            w:    std::ptr::null_mut(),
            vlen: Default::default(),
            dlen: Default::default(),
            elen: Default::default(),
            wlen: Default::default(),
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct SparseGraph {
    pub v: Vec<size_t>,
    pub d: Vec<::std::os::raw::c_int>,
    pub e: Vec<::std::os::raw::c_int>,
}

impl SparseGraph {
    pub fn new(vertices: usize, edges: usize) -> Self {
        SparseGraph {
            v: vec![0; vertices],
            d: vec![0; vertices],
            e: vec![0; edges],
        }
    }
}

impl<'a> std::convert::From<&'a mut SparseGraph> for sparsegraph {
    fn from(g: &'a mut SparseGraph) -> Self {
        sparsegraph {
            nv: g.v.len() as ::std::os::raw::c_int,
            nde: g.e.len() as size_t,
            v: g.v.as_mut_ptr(),
            d: g.d.as_mut_ptr(),
            e: g.e.as_mut_ptr(),
            w: std::ptr::null_mut(),
            vlen: g.v.len() as size_t,
            dlen: g.d.len() as size_t,
            elen: g.e.len() as size_t,
            wlen: 0,
        }
    }
}
