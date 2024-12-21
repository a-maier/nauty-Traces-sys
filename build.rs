use std::collections::HashMap;
use std::env;
use std::ffi::c_long;
use std::fmt::{Display, Write};
use std::path::PathBuf;

fn main() {
    let defines = get_nauty_defines();

    #[cfg(feature = "bundled")]
    compile_nauty(&defines);

    let wrapper = if cfg!(feature = "bundled") {
        println!("cargo:rustc-link-lib=nauty_bundled");
        "wrapper-bundled.h"
    } else {
        println!("cargo:rustc-link-lib=nauty");
        "wrapper.h"
    };
    println!("cargo:rerun-if-changed={wrapper}");

    let bindings = bindgen::Builder::default()
        .header_contents("nauty-flags.h", &c_defines(&defines))
        .header(wrapper)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .derive_eq(true)
        .derive_hash(true)
        .derive_ord(true)
        .no_partialeq("sparsegraph")
        .allowlist_type("optionblk")
        .allowlist_type("TracesOptions")
        .allowlist_type("TracesStats")
        .allowlist_type("statsblk")
        .allowlist_function("densenauty")
        .allowlist_function("nauty")
        .allowlist_function("Traces")
        .allowlist_function("nextelement")
        .allowlist_function("permset")
        .allowlist_function("orbjoin")
        .allowlist_function("writeperm")
        .allowlist_function("isautom")
        .allowlist_function("updatecan")
        .allowlist_function("refine")
        .allowlist_function("refine1")
        .allowlist_function("writegroupsize")
        .allowlist_function("nauty_check")
        .allowlist_function("groupautomproc")
        .allowlist_function("grouplevelproc")
        .allowlist_function("groupptr")
        .allowlist_function("makecosetreps")
        .allowlist_function("allgroup")
        .allowlist_function("sparsenauty")
        .allowlist_function("aresame_sg")
        .allowlist_function("sortlists_sg")
        .allowlist_function("put_sg")
        .allowlist_function("copy_sg")
        .allowlist_function("sg_to_nauty")
        .allowlist_function("nauty_to_sg")
        .allowlist_function("freeschreier")
        .allowlist_function("addpermutation")
        .allowlist_function("numcomponents")
        .allowlist_function("sources_sinks")
        .allowlist_function("digoncount")
        .allowlist_function("numind3sets1")
        .allowlist_function("numsquares")
        .allowlist_function("breakcellwt")
        .allowlist_function("settolist")
        .allowlist_function("listtoset")
        .allowlist_function("ran_init")
        .allowlist_function("ran_init_time")
        .allowlist_function("ran_nextran")
        .allowlist_function("adjacencies")
        .allowlist_function("adjacencies_sg")
        .allowlist_var("dispatch_graph")
        .allowlist_var("dispatch_sparse")
        // types are off for the following
        .allowlist_var("TRUE")
        .allowlist_var("FALSE")
        .allowlist_var("CONSOLWIDTH")
        .allowlist_var("MTOOBIG")
        .allowlist_var("NTOOBIG")
        .allowlist_var("CANONGNIL")
        .allowlist_var("NAUABORTED")
        .allowlist_var("NAUKILLED")
        // and this is completely off (linker error)
        // .allowlist_var("bit")
        .allowlist_var("NAUTYVERSIONID")
        .allowlist_var("WORDSIZE")
        .allowlist_var("stderr")
        .allowlist_var("stdout")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

#[cfg(feature = "bundled")]
fn compile_nauty<K: AsRef<str>, V: AsRef<str>>(
    defines: &HashMap<K, Option<V>>,
) {
    const NAUTY_DIR: &str = "nauty2_8_9";
    const NAUTY_HEADERS: [&str; 16] = [
        "gtools.h",
        "gutils.h",
        "naugroup.h",
        "naugstrings.h",
        "naurng.h",
        "nausparse.h",
        "nautaux.h",
        "nautinv.h",
        "naututil.h",
        "nautycliquer.h",
        "nauty.h",
        "planarity.h",
        "rng.h",
        "schreier.h",
        "traces.h",
        "sorttemplates.c",
    ];
    const NAUTY_SRC_FILES: [&str; 15] = [
        "nauty.c",
        "nautil.c",
        "nausparse.c",
        "naugraph.c",
        "schreier.c",
        "naurng.c",
        "traces.c",
        "gtools.c",
        "naututil.c",
        "nautinv.c",
        "gutil1.c",
        "gutil2.c",
        "gtnauty.c",
        "naugroup.c",
        "nautycliquer.c",
    ];

    for file in NAUTY_HEADERS.iter().chain(NAUTY_SRC_FILES.iter()) {
        println!("cargo:rerun-if-changed=src/{NAUTY_DIR}/{file}");
    }

    let mut cc_cmd = cc::Build::new();
    #[cfg(feature = "native")]
    cc_cmd.flag_if_supported("-march=native");
    #[cfg(feature = "popcnt")]
    {
        cc_cmd.flag_if_supported("-mpopcnt");
    }
    for (key, val) in defines.iter() {
        cc_cmd.define(key.as_ref(), val.as_ref().map(|v| v.as_ref()));
    }

    cc_cmd.warnings(false).files(
        NAUTY_SRC_FILES
            .iter()
            .map(|f| PathBuf::from_iter(["src", NAUTY_DIR, f])),
    );

    cc_cmd.compile("nauty_bundled");
}

fn get_nauty_defines() -> HashMap<&'static str, Option<&'static str>> {
    use std::env;

    if !cfg!(feature = "bundled") {
        return HashMap::new();
    }
    let mut defines = HashMap::new();

    const WORDSIZE_ENV_NAME: &str = "NAUTY_TRACES_WORDSIZE";
    const WORDSIZE_VAR_NAME: &str = "WORDSIZE";
    match env::var(WORDSIZE_ENV_NAME) {
        Ok(size) => {
            defines.insert(WORDSIZE_VAR_NAME, Some(&*size.leak()));
        }
        Err(env::VarError::NotPresent) => if c_long::BITS > 32 {
            // see section 3 in nauty manual (nauty 2.7r4)
            // TODO: this might be unnecessary, nauty sets it internally
            defines.insert(WORDSIZE_VAR_NAME, Some("64"));
        }
        Err(e) => panic!("Failed to fetch value of environment variable {WORDSIZE_ENV_NAME}: {e}")
    }

    const MAXN_ENV_NAME: &str = "NAUTY_TRACES_MAXN";
    const MAXN_VAR_NAME: &str = "MAXN";
    match env::var(MAXN_ENV_NAME) {
        Ok(size) => {
            defines.insert(MAXN_VAR_NAME, Some(&*size.leak()));
        }
        Err(env::VarError::NotPresent) => { }
        Err(e) => panic!("Failed to fetch value of environment variable {MAXN_ENV_NAME}: {e}")
    }

    #[cfg(feature = "tls")]
    defines.insert("USE_TLS", None);

    #[cfg(feature = "popcnt")]
    if is_x86_feature_detected!("popcnt") {
        defines.insert("HAVE_HWPOPCNT", Some("1"));
    }

    #[cfg(feature = "lzc")]
    if is_x86_feature_detected!("lzcnt") {
        defines.insert("HAVE_HWLZCNT", Some("1"));
    }

    defines
}

fn c_defines<K: Display, V: Display>(
    defines: &HashMap<K, Option<V>>,
) -> String {
    let mut define_str = String::new();
    for (key, val) in defines.iter() {
        if let Some(val) = val {
            writeln!(&mut define_str, "#define {key} {val}").unwrap();
        } else {
            writeln!(&mut define_str, "#define {key}").unwrap();
        }
    }
    define_str
}
