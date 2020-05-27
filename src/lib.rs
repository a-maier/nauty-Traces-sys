mod generated_bindings;
pub mod bindings;

#[cfg(test)]
mod tests {

    use super::bindings::*;
    use ::std::os::raw::c_int;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    // test nauty examples nautyex2.c to nautyex10.c

    #[test]
    fn nautyex2() {
        log_init();

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
        log_init();

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
    #[ignore] // fails for unknown reason
    fn nautyex4() {
        log_init();

        let n_range = 1..20;

        use ::std::os::raw::c_int;
        let mut options = optionblk::default();
        let mut stats = statsblk::default();
        options.writeautoms = TRUE;

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

            print!("Automorphism group size = ");
            unsafe {
                writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
            }
            println!();
        }
    }

    #[test]
    #[ignore] // same problem as nautyex4
    fn nautyex5() {
        log_init();

        let n_range = (2..20).step_by(2);

        let mut options = optionblk::default();
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

            if unsafe{ aresame_sg(&mut cg1,&mut cg2) } == TRUE {
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
        log_init();

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
        log_init();

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
        log_init();

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
        log_init();

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
        log_init();

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
