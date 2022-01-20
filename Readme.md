# Low-level bindings for nauty and Traces

This crate provides ffi bindings for
[nauty and Traces](http://pallini.di.uniroma1.it/).
Nauty and Traces are implementations of algorithms for computing graph
automorphisms.

# Usage

Add the following lines to your Cargo.toml:

```toml
[dependencies]
nauty-Traces-sys = "0.2"
```

# Caveats

* Nauty and Traces are not included in the sources and have to be
  already installed on your system.

* This crate was only tested against versions 2.6 and 2.7 of nauty and
  Traces.

* Some C macros have no direct equivalent.

    - Instead of using `DYNALLSTAT` and `DYNALLOC` you can create
      `Vec`s or arrays.

    - Use the `empty_graph` function in place of `EMPTY_GRAPH`.

    - Most `DEFAULT`-type macros have been replaced with
      implementations of the rust `Default` trait. For `sparsenauty`
      default options use `optionsblk::default_sparse()`.

    - The `SparseGraph` struct helps with creating sparse graphs. A
      `&mut SparseGraph` can be converted to the `sparsegraph` used by
      nauty and Traces.

* Sometimes you may want to free memory that was allocated internally
  by nauty and Traces, for example by the `nauty_to_sg` function. This
  can only be done safely if nauty and Traces are linked to the same
  libc as this crate. If that is the case, you can enable bindings to
  `DYNFREE` and `SG_FREE` by adding the following lines to your
  Cargo.toml:

  ```toml
  [dependencies]
  nauty-Traces-sys = { version = "0.2", features = ["libc"] }
  ```

# Examples

The following program prints the generators for the automorphism
groups of n-vertex polygons. It is a pretty literal translation of the
`nautyex2` C program that is part of the nauty and Traces bundle.

```rust
use nauty_Traces_sys::*;
use std::io::{self, Write};
use std::os::raw::c_int;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut options = optionblk::default();
    options.writeautoms = TRUE;
    let mut stats = statsblk::default();

    loop {
        print!("\nenter n : ");
        io::stdout().flush().unwrap();
        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        let n = input.trim().parse()?;
        if n > 0 {

            let m = SETWORDSNEEDED(n);

            unsafe {
                nauty_check(WORDSIZE as c_int, m as c_int, n as c_int, NAUTYVERSIONID as c_int);
            }

            let mut lab = vec![0; n];
            let mut ptn = vec![0; n];
            let mut orbits = vec![0; n];

            let mut g = empty_graph(m, n);
            for v in 0..n {
                ADDONEEDGE(&mut g, v, (v + 1) % n, m);
            }

            println!("Generators for Aut(C[{}]):", n);

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

            print!("[");
            for orbit in orbits {
                print!("{} ", orbit)
            }
            println!("]");

            print!("order = ");
            io::stdout().flush().unwrap();
            unsafe {
                writegroupsize(stderr, stats.grpsize1, stats.grpsize2);
            }
            println!();
        } else {
            break;
        }
    }
    Ok(())
}
```
