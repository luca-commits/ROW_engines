// Enable OpenCASCADE geometry kernel
SetFactory("OpenCASCADE");
// Define scale factor
scale = 0.01; // Scaling down by factor of 100
// Define the single outer rectangle (Transformer Core)
outer_rect = newv;
Rectangle(outer_rect) = {-30*scale, -30*scale, 0, 60*scale, 60*scale};
// Define the two medium rectangles (Air Region between coils and core)
left_medium_rect = newv;
Rectangle(left_medium_rect) = {-25*scale, -20*scale, 0, 20*scale, 40*scale};
right_medium_rect = newv;
Rectangle(right_medium_rect) = {5*scale, -20*scale, 0, 20*scale, 40*scale};
// Define the two small rectangles inside the left medium rectangle (Coils)
left_small_rect1 = newv;
Rectangle(left_small_rect1) = {-23.5*scale, -15*scale, 0, 3.5*scale, 30*scale}; // Secondary coil (tag 6)
left_small_rect2 = newv;
Rectangle(left_small_rect2) = {-10.5*scale, -15*scale, 0, 3.5*scale, 30*scale}; // Primary coil (tag 2)
// Define the two small rectangles inside the right medium rectangle (Coils)
right_small_rect1 = newv;
Rectangle(right_small_rect1) = {7*scale, -15*scale, 0, 3.5*scale, 30*scale}; // Primary coil (tag 5)
right_small_rect2 = newv;
Rectangle(right_small_rect2) = {20*scale, -15*scale, 0, 3.5*scale, 30*scale}; // Secondary coil (tag 6)
// Define a large circle that encloses the outer rectangle
r_outer = 50 * scale;
disk_outer = newv;
Disk(disk_outer) = {0, 0, 0, r_outer, r_outer};
// Create the airbox (Outer Circle minus Transformer Core)
circle_region() = BooleanDifference{ Surface{disk_outer}; Delete; }{ Surface{outer_rect}; };
// Remove medium rectangles from the transformer core, making them part of air
core_region() = BooleanDifference{ Surface{outer_rect}; Delete; }{ Surface{left_medium_rect, right_medium_rect}; };
// Remove the small rectangles from their respective medium rectangles
left_air_region() = BooleanDifference{ Surface{left_medium_rect}; Delete; }{
Surface{left_small_rect1, left_small_rect2};
};
right_air_region() = BooleanDifference{ Surface{right_medium_rect}; Delete; }{
Surface{right_small_rect1, right_small_rect2};
};
// Assign physical surfaces with correct tags
Physical Surface(4) = {circle_region[], left_air_region[], right_air_region[]}; // Airbox (including the space between core and coils)
Physical Surface(3) = {core_region[]}; // Transformer Core
Physical Surface(6) = {left_small_rect1, right_small_rect2}; // Secondary coils
Physical Surface(2) = {left_small_rect2}; // Primary coil (outward current)
Physical Surface(5) = {right_small_rect1}; // Primary coil (inward current)

// Define points for mesh refinement at corners
Point(6) = {-5*scale, -20*scale, 0, 0.1*scale};   // Bottom right corner of left medium rectangle
Point(7) = {-5*scale, 20*scale, 0, 0.1*scale};    // Top right corner of left medium rectangle
Point(9) = {5*scale, -20*scale, 0, 0.1*scale};    // Bottom left corner of right medium rectangle
Point(12) = {5*scale, 20*scale, 0, 0.1*scale};    // Top left corner of right medium rectangle

// Field for mesh size control of secondary coils
Field[1] = Box;
Field[1].VIn = 0.5*scale;  // Very small mesh size inside secondary coils
Field[1].VOut = 4*scale;   // Larger mesh size outside
Field[1].XMin = -24*scale;
Field[1].XMax = -20*scale;
Field[1].YMin = -15*scale;
Field[1].YMax = 15*scale;

Field[2] = Box;
Field[2].VIn = 0.5*scale;
Field[2].VOut = 4*scale;
Field[2].XMin = 20*scale;
Field[2].XMax = 24*scale;
Field[2].YMin = -15*scale;
Field[2].YMax = 15*scale;

// Field for corner refinement
Field[3] = Distance;
Field[3].NodesList = {6, 7, 9, 12};

// Threshold field for corners
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 1*scale;    // Min mesh size at the corners
Field[4].LcMax = 1*scale;    // Max mesh size far from corners
Field[4].DistMin = 0.2*scale;   // Distance from corner where min mesh size starts
Field[4].DistMax = 5*scale;     // Distance from corner where max mesh size starts

// Combine the secondary coil box fields
Field[5] = Min;
Field[5].FieldsList = {1, 2};

// Combine all refinement fields
Field[6] = Min;
Field[6].FieldsList = {5, 4};

// Use the combined field to define the mesh size
Background Field = 6;

// General mesh settings
Mesh.CharacteristicLengthMin = 0.1*scale;
Mesh.CharacteristicLengthMax = 0.7*scale;
Mesh.Algorithm = 6;           // Frontal-Delaunay for 2D meshes
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Generate the mesh
//Mesh 2;
