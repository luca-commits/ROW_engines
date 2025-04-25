Include "parameters.geo";

Point(1) = {0, 0, 0, 1};

Point(6) = {0, radius_airgap_stator, 0, lc_airgap};
Point(7) = {radius_airgap_stator, 0, 0, lc_airgap};
Point(8) = {0, -radius_airgap_stator, 0, lc_airgap};
Point(9) = {-radius_airgap_stator, 0, 0, lc_airgap};

// Define points and curves for outer circle
Point(139) = {0, radius_airgap_rotor, 0, lc_airgap};
Point(140) = {radius_airgap_rotor, 0, 0, lc_airgap};
Point(141) = {0, -radius_airgap_rotor, 0, lc_airgap};
Point(142) = {-radius_airgap_rotor, 0, 0, lc_airgap};
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
