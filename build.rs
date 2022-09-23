use std::env;
use std::path::PathBuf;

fn main() {

    #[cfg(feature = "bundled")]
    compile_nauty();

    let wrapper;
    if cfg!(feature = "bundled") {
        println!("cargo:rustc-link-lib=nauty_bundled");
        wrapper = "wrapper-bundled.h";
    } else {
        println!("cargo:rustc-link-lib=nauty");
        wrapper = "wrapper.h";
    };
    println!("cargo:rerun-if-changed={wrapper}");

    let bindings = bindgen::Builder::default()
        .header(wrapper)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
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
fn compile_nauty() {
    const NAUTY_DIR: &str = "nauty27r4";
    const NAUTY_HEADERS: [&str; 5] = [
        "nauty-h.in",
        "schreier.h",
        "naugroup.h",
        "nausparse.h",
        "traces.h",
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
    cc_cmd.warnings(false).files(
        NAUTY_SRC_FILES.iter().map(|f| PathBuf::from_iter(["src", NAUTY_DIR, f]))
    );

    cc_cmd.compile("nauty_bundled");

}
