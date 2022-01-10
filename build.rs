use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-lib=nauty");
    println!("cargo:rerun-if-changed=wrapper.h");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
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
