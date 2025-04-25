// radius_cylinder = 1;
// airgap_stator_width = 0.1;
// airgap_rotor_width = 0.1; 
// radius_airgap_stator = radius_cylinder + airgap_stator_width;
// radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;

// lc_airgap = 0.1;

Include "parameters.geo";

Point(1) = {0, 0, 0, 1};

Point(6) = {0, radius_airgap_stator, 0, 0.05};
Point(7) = {radius_airgap_stator, 0, 0, 0.05};
Point(8) = {0, -radius_airgap_stator, 0, 0.05};
Point(9) = {-radius_airgap_stator, 0, 0, 0.05};

// Define points and curves for outer circle
Point(139) = {0, radius_airgap_rotor, 0, 0.05};
Point(140) = {radius_airgap_rotor, 0, 0, 0.05};
Point(141) = {0, -radius_airgap_rotor, 0, 0.05};
Point(142) = {-radius_airgap_rotor, 0, 0, 0.05};
Circle(1001) = {139, 1, 140};
//+
Circle(1002) = {140, 1, 141};
//+
Circle(1003) = {141, 1, 142};
//+
Circle(1004) = {142, 1, 139};

// Circles for outer boundary
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Curve Loop(2) = {5, 6, 7, 8};  // Outer loop
//+
Curve Loop(3) = {1004, 1001, 1002, 1003};
//+
Plane Surface(9) = {2, 3};
//+
Physical Surface("airgap", 4) = {9};
