cp ../../build/projects/cylinder_test/benchmark_file_rosenbrock-wanner_transformator_magnetic_field.csv .
cp ../../build/projects/cylinder_test/benchmark_file_rosenbrock-wanner_transformator_current.csv .
# cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_small_timestep_transformator_B_field.csv . 
# cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_small_timestep_transformator_B_power.csv . 
# cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_small_timestep_transformator_current.csv . 
cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_transformator_B_field.csv . 
cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_transformator_current.csv . 
cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_transformator_B_power.csv . 
cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_transformator_current.csv . 
cp ../../build/projects/cylinder_test/benchmark_file_rosenbrock-wanner_transformator_magnetic_power.csv . 
echo "\nL2 norm B field energy, BDF2 vs benchmark"
python3 compute_L2_norm.py benchmark_file_rosenbrock-wanner_transformator_magnetic_field.csv benchmark_file_newton_methods_small_timestep_transformator_B_field.csv
echo "\nL2 norm B field energy, ROW vs benchmark"
python3 compute_L2_norm.py benchmark_file_newton_methods_transformator_B_field.csv benchmark_file_newton_methods_small_timestep_transformator_B_field.csv
echo "\nL2 norm current, ROW vs Benchmark"
python3 compute_L2_norm.py benchmark_file_rosenbrock-wanner_transformator_current.csv benchmark_file_newton_methods_small_timestep_transformator_current.csv
echo "\nL2 norm current, BDF2 vs Benchmark"
python3 compute_L2_norm.py  benchmark_file_newton_methods_transformator_current.csv benchmark_file_newton_methods_small_timestep_transformator_current.csv
echo "\nL2 norm total power loss, ROW vs Benchmark"
python3 compute_L2_norm.py benchmark_file_rosenbrock-wanner_transformator_magnetic_power.csv benchmark_file_newton_methods_small_timestep_transformator_B_power.csv
echo "\nL2 norm total power loss, BDF2 vs Benchmark"
python3 compute_L2_norm.py  benchmark_file_newton_methods_transformator_B_power.csv benchmark_file_newton_methods_small_timestep_transformator_B_power.csv