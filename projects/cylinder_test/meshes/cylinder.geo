// Gmsh .geo file for 2D cross-section of a conducting cylinder in an airbox

// Parameters
radius_cylinder = 1;  // Radius of the conducting cylinder
length_airbox = 50;    // Length of the airbox
height_airbox = 50;    // Height of the airbox

// Define points (in the 2D plane)
Point(1) = {0, 0, 0, 1.0};                        // Center of the conducting cylinder
Point(2) = {length_airbox/2, height_airbox/2, 0, 2};  // Top-right corner of the airbox
Point(3) = {-length_airbox/2, height_airbox/2, 0, 2}; // Top-left corner of the airbox
Point(4) = {-length_airbox/2, -height_airbox/2, 0, 2}; // Bottom-left corner of the airbox
Point(5) = {length_airbox/2, -height_airbox/2, 0, 2};  // Bottom-right corner of the airbox

// Define the lines for the airbox
Line(7) = {2, 3};  // Top of the airbox
Line(8) = {3, 4};  // Left side of the airbox
Line(9) = {4, 5};  // Bottom of the airbox
Line(10) = {5, 2}; // Right side of the airbox

Curve Loop(1) = {7, 8 , 9, 10};

// Define surfaces
//Plane Surface(12) = {1}; // Surface of the airbox

// Center of the conducting cylinder
Point(11) = {0, radius_cylinder, 0, 0.05};  // Top-right corner of the airbox //0.02
Point(12) = {radius_cylinder, 0, 0, 0.05};  // Top-left corner of the airbox
Point(13) = {0, -radius_cylinder, 0, 0.05}; // Bottom-left corner of the airbox
Point(14) = {-radius_cylinder, 0, 0, 0.05}; // Bottom-right corner of the airbox

Circle(15) = {11, 1, 12};
Circle(16) = {12, 1, 13};
Circle(17) = {13, 1, 14};
Circle(18) = {14, 1, 11};

Line Loop(2) = {15, 16, 17, 18};

Plane Surface(2) = {1, 2};
Plane Surface(3) = {2};


Physical Surface(1) = {2};
Physical Surface(2) = {3};

// Mesh generation
Mesh.MeshSizeMax = 100;//0.3;
Mesh 2;
// Mesh.SaveElementTagType = 2;
