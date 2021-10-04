#!/bin/sh

RUST_LOG=info target/debug/gis-algorithms -e 1 -ncp 10
fish diagrams.fish
