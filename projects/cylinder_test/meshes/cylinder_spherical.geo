// Gmsh .geo file for 2D cross-section of a conducting cylinder in an airbox

// Parameters
radius_cylinder = 1;  // Radius of the conducting cylinder
radius_airbox = 20; 
// Define points (in the 2D plane)
lc_airbox = 0.2;
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

lc = 0.03;
// Center of the conducting cylinder
Point(11) = {0, radius_cylinder, 0, lc};  
Point(12) = {radius_cylinder, 0, 0, lc};  
Point(13) = {0, -radius_cylinder, 0, lc}; 
Point(14) = {-radius_cylinder, 0, 0, lc}; 

Circle(15) = {11, 1, 12};
Circle(16) = {12, 1, 13};
Circle(17) = {13, 1, 14};
Circle(18) = {14, 1, 11};

Line Loop(2) = {15, 16, 17, 18};

Plane Surface(2) = {1, 2};
Plane Surface(3) = {2};


Physical Surface(1) = {2};
Physical Surface(2) = {3};

Mesh 2;


