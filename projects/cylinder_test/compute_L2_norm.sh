cp ../../build/projects/cylinder_test/benchmark_file_newton_methods_box_irregular_ring.csv .
cp ../../build/projects/cylinder_test/benchmark_file_rosenbrock-wanner_box_irregular_ring.csv .
python3 compute_L2_norm.py benchmark_file_rosenbrock-wanner_box_irregular_ring.csv benchmark_file_newton_methods_box_irregular_ring.csv 