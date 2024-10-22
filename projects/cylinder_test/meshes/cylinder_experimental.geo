// Parameters
radius_inner_cylinder = 0.4;         // Radius of the inner cylinder
radius_outer_conducting = 1;       // Radius of the outer conducting cylinder
radius_outer_air = 1.1;              // Radius of the outer air cylinder
length_airbox = 5;                 // Half-length of the square airbox
height_airbox = 5;                 // Half-height of the square airbox

// Mesh sizes for different regions
lc_inner_cylindre_conducting = 0.06 ;                      // Coarse mesh for inner cylinder
lc_cylindre_conducting = 0.02;         // Fine mesh for outer conducting cylinder
lc_outer_cylinder_air = 0.08;                // Fine mesh for outer air cylinder
lc_airbox = 1;

// Define the center point
Point(0) = {0, 0, 0, 1};           // Center point for all circles

// Define the inner cylinder
Point(5) = {0, radius_inner_cylinder, 0, lc_inner_cylindre_conducting};
Point(6) = {radius_inner_cylinder, 0, 0, lc_inner_cylindre_conducting};
Point(7) = {0, -radius_inner_cylinder, 0, lc_inner_cylindre_conducting};
Point(8) = {-radius_inner_cylinder, 0, 0, lc_inner_cylindre_conducting};

Circle(9) = {5, 0, 6};
Circle(10) = {6, 0, 7};
Circle(11) = {7, 0, 8};
Circle(12) = {8, 0, 5};

Line Loop(1) = {9, 10, 11, 12};

// Define the outer conducting cylinder
Point(13) = {0, radius_outer_conducting, 0, lc_cylindre_conducting};
Point(14) = {radius_outer_conducting, 0, 0, lc_cylindre_conducting};
Point(15) = {0, -radius_outer_conducting, 0, lc_cylindre_conducting};
Point(16) = {-radius_outer_conducting, 0, 0, lc_cylindre_conducting};

Circle(17) = {13, 0, 14};
Circle(18) = {14, 0, 15};
Circle(19) = {15, 0, 16};
Circle(20) = {16, 0, 13};

Line Loop(2) = {17, 18, 19, 20};

// Define the outer air cylinder
Point(21) = {0, radius_outer_air, 0, lc_outer_cylinder_air};
Point(22) = {radius_outer_air, 0, 0, lc_outer_cylinder_air};
Point(23) = {0, -radius_outer_air, 0, lc_outer_cylinder_air};
Point(24) = {-radius_outer_air, 0, 0, lc_outer_cylinder_air};

Circle(25) = {21, 0, 22};
Circle(26) = {22, 0, 23};
Circle(27) = {23, 0, 24};
Circle(28) = {24, 0, 21};

Line Loop(3) = {25, 26, 27, 28};

// Define the airbox boundary points
Point(1) = {length_airbox, height_airbox, 0, lc_airbox};
Point(2) = {-length_airbox, height_airbox, 0, lc_airbox};
Point(3) = {-length_airbox, -height_airbox, 0, lc_airbox};
Point(4) = {length_airbox, -height_airbox, 0, lc_airbox};

Line(29) = {1, 2};
Line(30) = {2, 3};
Line(31) = {3, 4};
Line(32) = {4, 1};

Line Loop(4) = {29, 30, 31, 32};

// Define surfaces
Plane Surface(1) = {4, 3};
Plane Surface(2) = {3, 2};
Plane Surface(3) = {2, 1};
Plane Surface(4) = {1};

// Define physical groups
Physical Surface(1) = {1, 2};
Physical Surface(2) = {3, 4};

Mesh 2;
