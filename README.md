# openbench-seq

The sequence engine of [OpenBench](https://github.com/aarteixeira/openbench): motif search with substitutions
(IUPAC, both strands, circular), restriction cut sites with Dam/Dcm/CpG overlap flags, ORFs across the origin,
translation with NCBI genetic codes, and SnapGene `.dna` and ABIF `.ab1` parsers. No dependencies; compiled to
WebAssembly and loaded by OpenBench's browser workers and its server.

```sh
cargo test --locked
cargo build --locked --release --target wasm32-unknown-unknown   # target/wasm32-unknown-unknown/release/openbench_seq.wasm
```

The toolchain is pinned in `rust-toolchain.toml`; OpenBench vendors the built `.wasm` and its CI rebuilds it from
the pinned commit and compares the bytes. The calling convention is described at the top of `src/lib.rs`.

Apache-2.0.
