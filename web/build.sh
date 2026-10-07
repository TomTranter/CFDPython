#!/usr/bin/env sh
# Build the Rust solver to WebAssembly and copy it next to index.html.
# Needs: rustup target add wasm32-unknown-unknown
set -e
cd "$(dirname "$0")"
cargo build --manifest-path cavity/Cargo.toml --release --target wasm32-unknown-unknown
cp cavity/target/wasm32-unknown-unknown/release/cavity.wasm cavity.wasm
echo "wrote web/cavity.wasm ($(wc -c < cavity.wasm) bytes)"
