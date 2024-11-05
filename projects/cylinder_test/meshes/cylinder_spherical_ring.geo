// Gmsh .geo file for 2D cross-section of a conducting cylinder in an airbox

// Parameters
radius_cylinder = 1;  // Radius of the conducting cylinder
radius_airbox = 10; 
inner_rad_ring = 1.4; 
outer_rad_ring = 1.45; 


lc_airbox = 3;
// Define points (in the 2D plane)
Point(1) = {0, 0, 0, 1};                        // Center of the conducting cylinder
Point(2) = {0, radius_airbox, 0, lc_airbox};  // Top-right corner of the airbox
Point(3) = {radius_airbox, 0, 0, lc_airbox}; // Top-left corner of the airbox
Point(4) = {0, -radius_airbox, 0, lc_airbox}; // Bottom-left corner of the airbox
Point(5) = {-radius_airbox, 0, 0, lc_airbox};  // Bottom-right corner of the airbox

Circle(7) = {2, 1, 3}; 
Circle(8) = {3, 1, 4};
Circle(9) = {4, 1, 5};
Circle(10) = {5, 1, 2};

Curve Loop(1) = {7, 8 , 9, 10};

// Define surfaces
//Plane Surface(12) = {1}; // Surface of the airbox

lc_cylinder = 0.1;
// Center of the conducting cylinder
Point(11) = {0, radius_cylinder, 0, lc_cylinder};  
Point(12) = {radius_cylinder, 0, 0, lc_cylinder};  
Point(13) = {0, -radius_cylinder, 0, lc_cylinder}; 
Point(14) = {-radius_cylinder, 0, 0, lc_cylinder}; 

Circle(15) = {11, 1, 12};
Circle(16) = {12, 1, 13};
Circle(17) = {13, 1, 14};
Circle(18) = {14, 1, 11};

Line Loop(2) = {15, 16, 17, 18}; //cylinder


lc_inner_ring = 0.1;
//inner circle ring
Point(19) = {0, inner_rad_ring, 0, lc_inner_ring};  
Point(20) = {inner_rad_ring, 0, 0, lc_inner_ring};  
Point(21) = {0, -inner_rad_ring, 0, lc_inner_ring}; 
Point(22) = {-inner_rad_ring, 0, 0, lc_inner_ring}; 

Circle(19) = {19, 1, 20};
Circle(20) = {20, 1, 21};
Circle(21) = {21, 1, 22};
Circle(22) = {22, 1, 19};

Line Loop(3) = {19, 20, 21, 22}; //inner ring circle


lc_outer_ring = 0.1;
//outer circle ring
Point(23) = {0, outer_rad_ring, 0, lc_outer_ring};  
Point(24) = {outer_rad_ring, 0, 0, lc_outer_ring};  
Point(25) = {0, -outer_rad_ring, 0, lc_outer_ring}; 
Point(26) = {-outer_rad_ring, 0, 0, lc_outer_ring}; 

Circle(23) = {23, 1, 24};
Circle(24) = {24, 1, 25};
Circle(25) = {25, 1, 26};
Circle(26) = {26, 1, 23};

Line Loop(4) = {23, 24, 25, 26}; //outer ring circle

Plane Surface(1) = {2};    //cylinder
Plane Surface(2) = {2, 3}; //cylinder-ring airgap
Plane Surface(3) = {3, 4}; //ring
Plane Surface(4) = {4, 1}; //airbox


Physical Surface(1) = {2, 4}; // tag 1: air
Physical Surface(2) = {1};    // inner cylinder
Physical Surface(3) = {3};    // ring

Mesh 2;
// Mesh.SaveElementTagType = 2;
