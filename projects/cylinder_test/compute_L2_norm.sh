#!/bin/zsh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 [transformator|box_irregular_ring]"
    exit 1
fi

GEOMETRY="$1"

if [ "$GEOMETRY" != "transformator" ] && [ "$GEOMETRY" != "box_irregular_ring" ]; then
    echo "Invalid geometry. Please use 'transformator' or 'box_irregular_ring'."
    exit 1
fi

# Define method suffixes
FILES=(
    "rosenbrock-wanner_${GEOMETRY}_magnetic_field"
    "rosenbrock-wanner_${GEOMETRY}_current"
    "rosenbrock-wanner_${GEOMETRY}_magnetic_power"
    "newton_methods_small_timestep_${GEOMETRY}_B_field"
    "newton_methods_small_timestep_${GEOMETRY}_current"
    "newton_methods_small_timestep_${GEOMETRY}_B_power"
    "newton_methods_${GEOMETRY}_B_field"
    "newton_methods_${GEOMETRY}_current"
    "newton_methods_${GEOMETRY}_B_power"
)

# Copy files
for f in "${FILES[@]}"; do
    cp "../../build/projects/cylinder_test/benchmark_file_${f}.csv" .
done

# Run comparisons
echo "\nL2 norm B field energy, ROW vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_rosenbrock-wanner_${GEOMETRY}_magnetic_field.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_B_field.csv"

echo "\nL2 norm B field energy, BDF2 vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_newton_methods_${GEOMETRY}_B_field.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_B_field.csv"

echo "\nL2 norm current, ROW vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_rosenbrock-wanner_${GEOMETRY}_current.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_current.csv"

echo "\nL2 norm current, BDF2 vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_newton_methods_${GEOMETRY}_current.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_current.csv"

echo "\nL2 norm total power loss, ROW vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_rosenbrock-wanner_${GEOMETRY}_magnetic_power.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_B_power.csv"

echo "\nL2 norm total power loss, BDF2 vs Benchmark"
python3 compute_L2_norm.py "benchmark_file_newton_methods_${GEOMETRY}_B_power.csv" "benchmark_file_newton_methods_small_timestep_${GEOMETRY}_B_power.csv"
