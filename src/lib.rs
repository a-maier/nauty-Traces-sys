#![doc(html_root_url = "https://docs.rs/nauty-Traces-sys/0.1.0")]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
mod bindings;

pub use bindings::*;
use num_integer::Integer;

// bindgen doesn't get the types right for the following constants
pub const FALSE: boolean = bindings::FALSE as boolean;
pub const TRUE: boolean = bindings::TRUE as boolean;
pub const CONSOLWIDTH: ::std::os::raw::c_int = bindings::CONSOLWIDTH as ::std::os::raw::c_int;

// bindgen gets this wrong somehow? linker won't find it.
#[allow(clippy::unreadable_literal)]
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

impl optionblk {
    pub fn default_sparse() -> Self {
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
            dispatch: unsafe {&mut dispatch_sparse},
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

/// Create an empty graph with `n` vertices.
/// `m` should be set to `SETWORDSNEEDED(n)`.
pub fn empty_graph(m: usize, n: usize) -> Vec<graph> {
    vec![0; m*n]
}

/// Create an uninitialised graph with `n` vertices.
/// `m` should be set to `SETWORDSNEEDED(n)`.
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

/// Sparse graph with allocated memory
#[derive(Debug, Default, Clone)]
pub struct SparseGraph {
    pub v: Vec<size_t>,
    pub d: Vec<::std::os::raw::c_int>,
    pub e: Vec<::std::os::raw::c_int>,
}

impl SparseGraph {
    /// Create a sparse graph with the given number of vertices and edges
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

#[cfg(test)]
mod tests {

    use super::*;
    use ::std::os::raw::c_int;

    // test nauty examples nautyex2.c to nautyex10.c

    #[test]
    fn nautyex2() {
        let orders = [1., 2., 6., 8., 10., 12., 14., 16., 18., 20.];

        let mut options = optionblk::default();
        // options.writeautoms = TRUE;
        let mut stats = statsblk::default();

        for (n, order) in orders.iter().enumerate() {
            let n = n + 1;
            let m = SETWORDSNEEDED(n);

            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab = Vec::with_capacity(n);
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut g = empty_graph(m, n);
            for v in 0..n {
                ADDONEEDGE(&mut g, v, (v + 1) % n, m);
            }
            unsafe {
                densenauty(
                    g.as_mut_ptr(),
                    lab.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    std::ptr::null_mut()
                );
            }
            assert_eq!(stats.grpsize1, *order);
            assert_eq!(stats.grpsize2, 0);
        }
    }

    extern "C" fn writeautom(p: *mut i32, n: i32) {
        for i in 0..n {
            print!(" {:2}", unsafe{ *p.offset(i as isize) })
        }
        println!()
    }

    #[test]
    fn nautyex3() {
        // nautyex3.c nauty example

        let mut options = optionblk::default();
        let mut stats = statsblk::default();

        /* The following cause nauty to call two procedures which
        store the group information as nauty runs. */

        options.userautomproc = Some(groupautomproc);
        options.userlevelproc = Some(grouplevelproc);

        for n in 1..20 {
            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int,n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut g = empty_graph(n, m);
            let mut lab = Vec::with_capacity(n);
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);

            for v in 0..n {
                ADDONEEDGE(&mut g, v, (v + 1) % n, m)
            }

            println!("Automorphisms of C[{}]:", n);

            unsafe {
                densenauty(
                    g.as_mut_ptr(),
                    lab.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    std::ptr::null_mut()
                );

                /* Get a pointer to the structure in which the group information
                has been stored.  If you use TRUE as an argument, the
                structure will be "cut loose" so that it won't be used
                again the next time nauty() is called.  Otherwise, as
                here, the same structure is used repeatedly. */

                let group = groupptr(FALSE);

                /* Expand the group structure to include a full set of coset
                representatives at every level.  This step is necessary
                if allgroup() is to be called. */

                makecosetreps(group);

                /* Call the procedure writeautom() for every element of the group.
                The first call is always for the identity. */

                allgroup(group, Some(writeautom));
            }
        }
    }

    #[test]
    fn nautyex4() {
        let n_range = 1..20;

        use ::std::os::raw::c_int;
        let mut options = optionblk::default_sparse();
        let mut stats = statsblk::default();
        // options.writeautoms = TRUE;

        for n in n_range {
            let m = SETWORDSNEEDED(n);

            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab = Vec::with_capacity(n);
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);

            /* SG_ALLOC makes sure that the v,d,e fields of a sparse graph
            structure point to arrays that are large enough.  This only
            works if the structure has been initialised. */

            let mut sg = SparseGraph::new(
                n,   /* Number of vertices */
                2*n  /* Number of directed edges */
            );

            for i in 0..n {
                sg.v[i] = 2*i as size_t;
                sg.d[i] = 2;
                sg.e[2*i] = ((i+n-1) % n) as c_int;      /* edge i->i-1 */
                sg.e[2*i+1] = ((i+n+1) % n) as c_int;    /* edge i->i+1 */
            }

            println!("Generators for Aut(C[{}]):", n);
            unsafe {
                sparsenauty(
                    &mut (&mut sg).into(),
                    lab.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    std::ptr::null_mut(),
                );
            }

            println!("Automorphism group size = {}", stats.grpsize1);
            if n < 3 {
                assert_eq!(stats.grpsize1, n as f64)
            } else {
                assert_eq!(stats.grpsize1, 2.*n as f64)
            }
            println!();
        }
    }

    #[test]
    fn nautyex5() {
        let n_range = (2..20).step_by(2);

        let mut options = optionblk::default_sparse();
        let mut stats = statsblk::default();

        let mut cg1 = sparsegraph::default();
        let mut cg2 = sparsegraph::default();

        /* Select option for canonical labelling */
        options.getcanon = TRUE;

        for n in n_range {
            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m  as c_int, n  as c_int, NAUTYVERSIONID  as c_int);
            }
            let mut lab1 = vec![0; n];
            let mut lab2 = vec![0; n];
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut map = vec![0; n];

            /* Now make the first graph */

            let mut sg1 = SparseGraph::new(
                n, /* Number of vertices */
                3*n /* Number of directed edges */
            );

            for i in 0..n {
                sg1.v[i] = (3*i) as size_t;     /* Position of vertex i in v array */
                sg1.d[i] = 3;                /* Degree of vertex i */
            }

            for i in (0..n).step_by(2) {
                sg1.e[sg1.v[i] as usize] = (i+1) as c_int;
                sg1.e[sg1.v[i+1] as usize] = i as c_int;
            }

            for i in 0..n-2 {  /* Clockwise edges */
                sg1.e[(sg1.v[i]+1) as usize] = (i+2) as c_int;
            }
            sg1.e[(sg1.v[n-2]+1) as usize] = 1;
            sg1.e[(sg1.v[n-1]+1) as usize] = 0;

            for i in 2..n  { /* Anticlockwise edges */
                sg1.e[(sg1.v[i]+2) as usize] = (i-2) as c_int;
            }
            sg1.e[(sg1.v[1]+2) as usize] = (n-2) as c_int;
            sg1.e[(sg1.v[0]+2) as usize] = (n-1) as c_int;

            /* Now make the second graph */

            let mut sg2 = SparseGraph::new(
                n,            /* Number of vertices */
                3*n,          /* Number of directed edges */
            );

            // // this is redundant already in nautyex5.c
            // for i in 0..n {
            //     sg2.v[i] = (3*i) as size_t;
            //     sg2.d[i] = 3;
            // }

            for i in 0..n {
                sg2.v[i] = (3*i) as size_t;
                sg2.d[i] = 3;
                sg2.e[sg2.v[i] as usize] = ((i+1) % n) as c_int;      /* Clockwise */
                sg2.e[(sg2.v[i]+1) as usize] = ((i+n-1) % n) as c_int;  /* Anti-clockwise */
                sg2.e[(sg2.v[i]+2) as usize] = ((i+n/2) % n) as c_int;  /* Diagonals */
            }

            /* Label sg1, result in cg1 and labelling in lab1; similarly sg2.
            It is not necessary to pre-allocate space in cg1 and cg2, but
            they have to be initialised as we did above.  */
            unsafe {
                sparsenauty(
                    &mut (&mut sg1).into(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg1,
                );
                sparsenauty(
                    &mut (&mut sg2).into(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg2,
                );
            }

            /* Compare canonically labelled graphs */
            let are_same = unsafe{ aresame_sg(&mut cg1,&mut cg2) } == TRUE;
            assert!(are_same);
            if are_same {
                println!("Isomorphic");
                if n <= 1000 {
                 /* Write the isomorphism.  For each i, vertex lab1[i]
                    of sg1 maps onto vertex lab2[i] of sg2.  We compute
                    the map in order of labelling because it looks better. */

                    for i in 0..n {
                        map[lab1[i] as usize] = lab2[i];
                    }
                    for i in 0..n {
                        print!(" {}-{}", i, map[i])
                    }
                    println!()
                }
            } else {
                println!("Not isomorphic.");
            }
        }
    }

    #[test]
    fn nautyex6() {
        let n_range = (2..20).step_by(2);

        let mut options = optionblk::default();
        let mut stats = statsblk::default();

        /* Select option for canonical labelling */

        options.getcanon = TRUE;

        for n in n_range {

            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab1 = vec![0; n];
            let mut lab2 = vec![0; n];
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut map = vec![0; n];
            let mut cg1 = empty_graph(n, m);
            let mut cg2 = empty_graph(n, m);

            /* Now make the first graph */
            /* ADDEDGE() is defined above */
            let mut g1 = empty_graph(m, n);

            for i in 0..n-2 {
                ADDONEEDGE(&mut g1, i, i+2, m)
            }
            ADDONEEDGE(&mut g1, n-2, 1, m);
            ADDONEEDGE(&mut g1, n-1, 0, m);
            for i in (0..n).step_by(2) {
                ADDONEEDGE(&mut g1,i,i+1,m)
            }

            /* Now make the second graph */

            let mut g2 = empty_graph(m, n);
            for i in 0..n-1 {
                ADDONEEDGE(&mut g2,i,i+1,m)
            }
            ADDONEEDGE(&mut g2,n-1,0,m);
            for i in 0..n/2 {
                ADDONEEDGE(&mut g2,i,i+n/2,m);
            }

            /* Label g1, result in cg1 and labelling in lab1; similarly g2.
            It is not necessary to pre-allocate space in cg1 and cg2, but
            they have to be initialised as we did above.  */

            unsafe {
                densenauty(
                    g1.as_mut_ptr(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    cg1.as_mut_ptr(),
                );

                densenauty(
                    g2.as_mut_ptr(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    cg2.as_mut_ptr(),
                );
            }

            /* Compare canonically labelled graphs */
            assert_eq!(cg1, cg2);
            if cg1 == cg2 {
                println!("Isomorphic.");
                if n <= 1000 {
                    /* Write the isomorphism.  For each i, vertex lab1[i]
                    of sg1 maps onto vertex lab2[i] of sg2.  We compute
                    the map in order of labelling because it looks better. */

                    for i in 0..n {
                        map[lab1[i] as usize] = lab2[i];
                    }
                    for i in 0..n {
                        print!(" {}-{}",i,map[i]);
                    }
                    println!()
                }
            }
        }
    }

    #[test]
    fn nautyex7() {
        let n_range = (2..20).step_by(2);

        let mut options = TracesOptions::default();
        let mut stats = TracesStats::default();

        /* Declare and initialize sparse graph structures */
        let mut cg1 = sparsegraph::default();
        let mut cg2 = sparsegraph::default();

        /* Select option for canonical labelling */

        options.getcanon = TRUE;

        for n in n_range {
            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab1 = vec![0; n];
            let mut lab2 = vec![0; n];
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut map = vec![0; n];

            /* Now make the first graph */
            let mut sg1 = SparseGraph::new(
                n,     /* Number of vertices */
                3*n    /* Number of directed edges */
            );

            for i in 0..n {
                sg1.v[i] = (3*i) as size_t;  /* Position of vertex i in v array */
                sg1.d[i] = 3;                /* Degree of vertex i */
            }

            for i in (0..n).step_by(2) { /* Spokes */
                sg1.e[sg1.v[i] as usize] = (i+1) as c_int;
                sg1.e[sg1.v[i+1] as usize] = i as c_int;
            }

            for i in 0..n-2 {  /* Clockwise edges */
                sg1.e[(sg1.v[i]+1) as usize] = (i+2) as c_int;
            }
            sg1.e[(sg1.v[n-2]+1) as usize] = 1;
            sg1.e[(sg1.v[n-1]+1) as usize] = 0;

            for i in 2..n {  /* Anticlockwise edges */
                sg1.e[(sg1.v[i]+2) as usize] = (i-2) as c_int;
            }
            sg1.e[(sg1.v[1]+2) as usize] = (n-2) as c_int;
            sg1.e[(sg1.v[0]+2) as usize] = (n-1) as c_int;

            /* Now make the second graph */
            let mut sg2 = SparseGraph::new(
                n,     /* Number of vertices */
                3*n    /* Number of directed edges */
            );

            for i in 0..n {
                sg2.v[i] = (3*i) as size_t;
                sg2.d[i] = 3;
                sg2.e[(sg2.v[i]) as usize] = ((i+1) % n) as c_int;      /* Clockwise */
                sg2.e[(sg2.v[i]+1) as usize] = ((i+n-1) % n) as c_int;  /* Anti-clockwise */
                sg2.e[(sg2.v[i]+2) as usize] = ((i+n/2) % n) as c_int;  /* Diagonals */
            }

         /* Label sg1, result in cg1 and labelling in lab1; similarly sg2.
            It is not necessary to pre-allocate space in cg1 and cg2, but
            they have to be initialised as we did above.  */

            unsafe {
                Traces(
                    &mut (&mut sg1).into(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg1
                );
                Traces(
                    &mut (&mut sg2).into(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg2
                );
            }


            /* Compare canonically labelled graphs */
            let are_same = unsafe{ aresame_sg(&mut cg1, &mut cg2) } == TRUE;
            assert!(are_same);
            if are_same {
                println!("Isomorphic.");
                if n <= 1000 {
                 /* Write the isomorphism.  For each i, vertex lab1[i]
                    of sg1 maps onto vertex lab2[i] of sg2.  We compute
                    the map in order of labelling because it looks better. */

                    for i in 0..n {
                        map[lab1[i] as usize] = lab2[i]
                    }
                    for i in 0..n {
                        print!(" {}-{}",i,map[i]);
                    }
                    println!();
                }
            }
            else {
                println!("Not isomorphic.");
            }
        }
    }

    #[test]
    fn nautyex8() {
        let n_range = (2..20).step_by(2);

        // nautyex8.c nauty example

        let mut options = optionblk::default();
        let mut stats = statsblk::default();

        options.getcanon = TRUE;

        for n in n_range {

            let m = SETWORDSNEEDED(n);

            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab1 = vec![0; n];
            let mut lab2 = vec![0; n];
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut map = vec![0; n];

            /* Now make the first graph */

            let mut g1 = empty_graph(m, n);

            for i in (0..n).step_by(2) { /* Spokes */
                ADDONEEDGE(&mut g1, i, i + 1, m);
            }

            for i in 0..n-2 { /* Cycle */
                ADDONEEDGE(&mut g1, i, i + 2, m);
            }
            ADDONEEDGE(&mut g1, 1, n - 2, m);
            ADDONEEDGE(&mut g1, 0, n - 1, m);

            /* Now make the second graph */

            let mut g2 = empty_graph(m, n);

            for i in 0..n {
                ADDONEEDGE(&mut g2, i, (i + 1) % n, m);     /* Rim */
                ADDONEEDGE(&mut g2, i, (i + n / 2) % n, m); /* Diagonals */
            }

            /* Create canonical graphs */

            let mut cg1 = empty_graph(m, n);
            let mut cg2 = empty_graph(m, n);

            unsafe {
                densenauty(
                    g1.as_mut_ptr(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    cg1.as_mut_ptr(),
                );

                densenauty(
                    g2.as_mut_ptr(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    m as c_int,
                    n as c_int,
                    cg2.as_mut_ptr(),
                );
            }

            assert_eq!(cg1, cg2);
            if cg1 == cg2 {
                println!("Isomorphic.");

                /* Write the isomorphism. For each i, vertex lab1[i]
                of sg1 maps onto vertex lab2[i] of sg2. We compute
                the map in order of labelling because it looks better. */
                for i in 0..n {
                    map[lab1[i] as usize] = lab2[i];
                }
                for (i, m) in map.iter().enumerate() {
                    print!(" {}-{}", i, m)
                }
                println!()
            } else {
                println!("Not isomorphic.")
            }
        }

    }

    #[test]
    fn nautyex9() {
        let mut auto_group_size = std::collections::HashMap::new();
        for i in &[3,4,6,7,8,9,11,12,14,15,16,18,19] {
            auto_group_size.insert(*i, None);
        }
        auto_group_size.insert(5, Some(10.));
        auto_group_size.insert(10, Some(320.));
        auto_group_size.insert(13, Some(78.));
        auto_group_size.insert(17, Some(136.));

        let n_range = 3..20;
        let mut options = TracesOptions::default();
        let mut stats = TracesStats::default();

        let mut gens = std::ptr::null_mut();

        /* Select option for passing generators to Traces */

        options.generators = &mut gens;

        for n in n_range {
            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab = Vec::with_capacity(n);
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut p = vec![0; n];
            let mut issquare = vec![FALSE; n];

            /* Initialise list of automorphisms */

            gens = std::ptr::null_mut();

            /* Find the squares and the degree */

            for i in 0..n {
                issquare[(i*i) % n] = TRUE;
            }
            if issquare.last() != Some(&TRUE) {
                assert_eq!(auto_group_size[&n], None);
                println!("-1 must be a square mod n; try again");
                continue;
            }

            let mut deg = 0;
            for i in 1..n {
                if issquare[i] == TRUE { deg += 1 }
            }

            /* Now make the graph */
            let mut sg = SparseGraph::new(
                n,     /* Number of vertices */
                n*deg    /* Number of directed edges */
            );

            for i in 0..n {
                sg.v[i] = (i*deg) as size_t;  /* Position of vertex i in v array */
                sg.d[i] = deg as c_int;      /* Degree of vertex i */
            }

            for i in 0..n {
                let mut k = sg.v[i] as usize;
                for j in 1..n {
                    if issquare[j] == TRUE {
                        sg.e[k] = ((i + j) % n) as c_int;
                        k += 1;
                    }
                }
            }

            /* Add known automorphisms */

            /* We wouldn't need freeschreier() if we were only
               processing one graph, but it doesn't hurt.  This
               is how to properly dispose of previous generators. */

            unsafe{ freeschreier(std::ptr::null_mut(), &mut gens); }

            /* Cyclic rotation */
            for i in 0..n {
                p[i] = ((i + 1) % n) as c_int;
            }
            unsafe{ addpermutation(&mut gens, p.as_mut_ptr(), n as c_int); }

            /* Reflection about 0 */
            for i in 0..n {
                p[i] = ((n - i) % n) as c_int;
            }
            unsafe{ addpermutation(&mut gens, p.as_mut_ptr(), n as c_int); }

             /* Call Traces */
            unsafe {
                Traces(
                    &mut (&mut sg).into(),
                    lab.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    std::ptr::null_mut(),
                );
            }

            assert_eq!(auto_group_size[&n], Some(stats.grpsize1));
            print!("Automorphism group size = ");
            unsafe {
                writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
            }
            println!();

            /* Traces left the automorphims we gave it, augmented by
            any extra automorphims it found, in a circular list
            pointed to by gens.  See schreier.txt for documentation. */
        }
    }

    #[test]
    fn nautyex10() {
        let n_range = (2..20).step_by(2);

        let mut options = TracesOptions::default();
        let mut stats = TracesStats::default();

        /* Declare and initialize sparse graph structures */
        let mut cg1 = sparsegraph::default();
        let mut cg2 = sparsegraph::default();


        for n in n_range {
            let m = SETWORDSNEEDED(n);
            unsafe {
                nauty_check(WORDSIZE as c_int, m  as c_int, n  as c_int, NAUTYVERSIONID  as c_int);
            }
            let mut lab1 = vec![0; n];
            let mut lab2 = vec![0; n];
            let mut ptn = Vec::with_capacity(n);
            let mut orbits = Vec::with_capacity(n);
            let mut map = vec![0; n];

            /* Now make the first graph */

            let mut sg1 = SparseGraph::new(
                n, /* Number of vertices */
                3*n /* Number of directed edges */
            );

            for i in 0..n {
                sg1.v[i] = (3*i) as size_t;     /* Position of vertex i in v array */
                sg1.d[i] = 3;                /* Degree of vertex i */
            }

            for i in (0..n).step_by(2) {
                sg1.e[sg1.v[i] as usize] = (i+1) as c_int;
                sg1.e[sg1.v[i+1] as usize] = i as c_int;
            }

            for i in 0..n-2 {  /* Clockwise edges */
                sg1.e[(sg1.v[i]+1) as usize] = (i+2) as c_int;
            }
            sg1.e[(sg1.v[n-2]+1) as usize] = 1;
            sg1.e[(sg1.v[n-1]+1) as usize] = 0;

            for i in 2..n  { /* Anticlockwise edges */
                sg1.e[(sg1.v[i]+2) as usize] = (i-2) as c_int;
            }
            sg1.e[(sg1.v[1]+2) as usize] = (n-2) as c_int;
            sg1.e[(sg1.v[0]+2) as usize] = (n-1) as c_int;

            /* Now make the second graph */

            let mut sg2 = SparseGraph::new(
                n,            /* Number of vertices */
                3*n,          /* Number of directed edges */
            );

            for i in 0..n {
                sg2.v[i] = (3*i) as size_t;
                sg2.d[i] = 3;
                sg2.e[sg2.v[i] as usize] = ((i+1) % n) as c_int;      /* Clockwise */
                sg2.e[(sg2.v[i]+1) as usize] = ((i+n-1) % n) as c_int;  /* Anti-clockwise */
                sg2.e[(sg2.v[i]+2) as usize] = ((i+n/2) % n) as c_int;  /* Diagonals */
            }

            /* Now we make the canonically labelled graphs by a two-step
            process.  The first call to Traces computes the
            automorphism group.  The second call computes the
            canonical labelling, using the automorphism group from
            the first call.

            We have declared a variable "generators" that will be
            used to hold the group generators between the two calls.
            It has to be initialised to NULL and its address has to
            be given to Traces using options.generators.  After the
            second call, we need to discard the generators with a
            call to freeschreier(), which also initializes it again. */

            let mut generators = std::ptr::null_mut();
            options.generators = &mut generators;

            options.getcanon = FALSE;
            unsafe {
                Traces(
                    &mut (&mut sg1).into(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    std::ptr::null_mut(),
                );
            }
            options.getcanon = TRUE;
            unsafe {
                Traces(
                    &mut (&mut sg1).into(),
                    lab1.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg1
                );
                freeschreier(std::ptr::null_mut(), &mut generators);
            }

            options.getcanon = FALSE;
            unsafe {
                Traces(
                    &mut (&mut sg2).into(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    std::ptr::null_mut(),
                );
            }
            options.getcanon = TRUE;
            unsafe {
                Traces(
                    &mut (&mut sg2).into(),
                    lab2.as_mut_ptr(),
                    ptn.as_mut_ptr(),
                    orbits.as_mut_ptr(),
                    &mut options,
                    &mut stats,
                    &mut cg2
                );
                freeschreier(std::ptr::null_mut(), &mut generators);
            }


            /* Compare canonically labelled graphs */
            let are_same = unsafe{ aresame_sg(&mut cg1,&mut cg2) } == TRUE;
            assert!(are_same);
            if are_same {
                println!("Isomorphic");
                if n <= 1000 {
                 /* Write the isomorphism.  For each i, vertex lab1[i]
                    of sg1 maps onto vertex lab2[i] of sg2.  We compute
                    the map in order of labelling because it looks better. */

                    for i in 0..n {
                        map[lab1[i] as usize] = lab2[i];
                    }
                    for i in 0..n {
                        print!(" {}-{}", i, map[i])
                    }
                    println!()
                }
            } else {
                println!("Not isomorphic.");
            }
        }
    }

}
