radius_airgap = 1.3;  // Radius of the conducting cylinder
lc_airbox = 0.1;
radius_airbox = 3;
radius_cylinder = 1;
airgap_stator_width = 0.1;
airgap_rotor_width = 0.1; 
radius_airgap_stator = radius_cylinder + airgap_stator_width;
radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;
lc_airgap = 0.05;
lc_cylinder = 0.1;
airgap_width = 0.2;
radius_ring_internal = radius_airgap + airgap_width;
ring_width = 0.3;
radius_ring_external = radius_ring_internal + ring_width;

// Mesh characteristic lengths
lc = 0.1;

// Center point
Point(1) = {0, 0, 0, lc};

// Points for air gap boundary (innermost)
Point(2) = {0, radius_airgap_rotor, 0, lc_airgap};
Point(3) = {radius_airgap_rotor, 0, 0, lc_airgap};
Point(4) = {0, -radius_airgap_rotor, 0, lc_airgap};
Point(5) = {-radius_airgap_rotor, 0, 0, lc_airgap};

// Points for ring inner boundary
Point(6) = {0, radius_ring_internal, 0, lc};
Point(7) = {radius_ring_internal, 0, 0, lc};
Point(8) = {0, -radius_ring_internal, 0, lc};
Point(9) = {-radius_ring_internal, 0, 0, lc};

// Points for ring outer boundary
Point(10) = {0, radius_ring_external, 0, lc};
Point(11) = {radius_ring_external, 0, 0, lc};
Point(12) = {0, -radius_ring_external, 0, lc};
Point(13) = {-radius_ring_external, 0, 0, lc};

// Points for outer air boundary
Point(14) = {0, radius_airbox, 0, lc};
Point(15) = {radius_airbox, 0, 0, lc};
Point(16) = {0, -radius_airbox, 0, lc};
Point(17) = {-radius_airbox, 0, 0, lc};

// Create circles
Circle(1) = {2, 1, 3};  // Innermost circle (airgap)
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {6, 1, 7};  // Ring inner boundary
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Circle(9) = {10, 1, 11};  // Ring outer boundary
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};

Circle(13) = {14, 1, 15};  // Outer air boundary
Circle(14) = {15, 1, 16};
Circle(15) = {16, 1, 17};
Circle(16) = {17, 1, 14};

// Create curve loops
Curve Loop(1) = {1, 2, 3, 4};    // Airgap boundary
Curve Loop(2) = {5, 6, 7, 8};    // Ring inner boundary
Curve Loop(3) = {9, 10, 11, 12}; // Ring outer boundary
Curve Loop(4) = {13, 14, 15, 16}; // Air outer boundary

// Create surfaces
Plane Surface(3) = {2, 1};    // Airgap region
Plane Surface(4) = {3, 2};    // Ring region
Plane Surface(5) = {4, 3};    // Outer air region


// Generate 2D mesh
Mesh 2;//+
Physical Surface("air", 1) = {3, 5};
//+
Physical Surface("ring", 3) = {4};


