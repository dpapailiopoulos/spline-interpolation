#!/bin/bash
set -e

progname="${0##*/}"
usage="Usage: $progname <spline_type>
Spline types:
  cubic     Run cubic spline interpolation
  bspline   Run B-spline interpolation
Options:
  -h        Show this help message
"

# Parse options
while getopts ':h' opt; do
  case "$opt" in
    h)  echo "$usage"; exit 0 ;;
    ?)  echo "$usage"; exit 1 ;;
  esac
done

# Get script and build directories
CURRENT_DIR=$(pwd)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
build_dir="$SCRIPT_DIR/build"

# Create build directory if missing
mkdir -p "$build_dir"

# Compile project
cmake -S "$SCRIPT_DIR" -B "$build_dir"
cmake --build "$build_dir"
cd "$build_dir"

# Determine interpolation type
case "$1" in
  cubic)
    output="cubic_spline_points.txt"
    "$build_dir/cubic_spline"
    ;;
  bspline)
    output="b_spline_points.txt"
    "$build_dir/b_spline"
    ;;
  *)
    echo "Unknown spline type: '$1'"
    echo "Expected 'cubic' or 'bspline'"
    exit 1
    ;;
esac

# Visualize results
python3 "$SCRIPT_DIR/visualize.py" "$build_dir/$output"

# Return to original dir
cd "$CURRENT_DIR"
