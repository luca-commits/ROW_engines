// Parameters
radius_cylinder = 1;  // Radius of the conducting cylinder
airgap_stator_width = 0.1;
airgap_rotor_width = 0.1; 
radius_airgap_stator = radius_cylinder + airgap_stator_width;
radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;

lc_airbox = 0.01;

Point(1) = {0, 0, 0, 1};
Point(2) = {0, radius_cylinder, 0, 0.1};
Point(3) = {radius_cylinder, 0, 0, 0.1};
Point(4) = {0, -radius_cylinder, 0, 0.1};
Point(5) = {-radius_cylinder, 0, 0, 0.1};

Point(6) = {0, radius_airgap_stator, 0, 0.05};
Point(7) = {radius_airgap_stator, 0, 0, 0.05};
Point(8) = {0, -radius_airgap_stator, 0, 0.05};
Point(9) = {-radius_airgap_stator, 0, 0, 0.05};

// Circles for inner boundary
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Circles for outer boundary
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Define curve loops
Curve Loop(1) = {1, 2, 3, 4};  // Inner loop
Curve Loop(2) = {5, 6, 7, 8};  // Outer loop

// Define surfaces
Plane Surface(1) = {1};        // Inner circle
Plane Surface(2) = {2, 1};     // Annular region

// Generate mesh

//+
Physical Surface("air", 1) = {2};
Physical Surface("cylinder", 2) = {1};
Mesh 2;//+
