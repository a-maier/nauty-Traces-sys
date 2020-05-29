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
        .whitelist_type("optionblk")
        .whitelist_type("TracesOptions")
        .whitelist_type("TracesStats")
        .whitelist_type("statsblk")
        .whitelist_function("densenauty")
        .whitelist_function("nauty")
        .whitelist_function("Traces")
        .whitelist_function("nextelement")
        .whitelist_function("permset")
        .whitelist_function("orbjoin")
        .whitelist_function("writeperm")
        .whitelist_function("isautom")
        .whitelist_function("updatecan")
        .whitelist_function("refine")
        .whitelist_function("refine1")
        .whitelist_function("writegroupsize")
        .whitelist_function("nauty_check")
        .whitelist_function("groupautomproc")
        .whitelist_function("grouplevelproc")
        .whitelist_function("groupptr")
        .whitelist_function("makecosetreps")
        .whitelist_function("allgroup")
        .whitelist_function("sparsenauty")
        .whitelist_function("aresame_sg")
        .whitelist_function("freeschreier")
        .whitelist_function("addpermutation")
        .whitelist_var("dispatch_graph")
        .whitelist_var("dispatch_sparse")
    // types are off for the following
        .whitelist_var("TRUE")
        .whitelist_var("FALSE")
        .whitelist_var("CONSOLWIDTH")
        .whitelist_var("MTOOBIG")
        .whitelist_var("NTOOBIG")
        .whitelist_var("CANONGNIL")
        .whitelist_var("NAUABORTED")
        .whitelist_var("NAUKILLED")
    // and this is completely off (linker error)
    // .whitelist_var("bit")
        .whitelist_var("NAUTYVERSIONID")
        .whitelist_var("WORDSIZE")
        .whitelist_var("stderr")
        .whitelist_var("stdout")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
