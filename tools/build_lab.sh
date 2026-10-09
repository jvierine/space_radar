#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "$0")/.." && pwd)
cargo build --release --target wasm32-unknown-unknown --manifest-path "$repo_root/lab-core/Cargo.toml"
cp "$repo_root/lab-core/target/wasm32-unknown-unknown/release/fmcw_lab_core.wasm" "$repo_root/web/lab/core.wasm"
chmod 644 "$repo_root/web/lab/core.wasm"
