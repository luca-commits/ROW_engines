// Parameters
// radius_cylinder = 1;  // Radius of the conducting cylinder
// airgap_stator_width = 0.1;
// airgap_rotor_width = 0.1; 
// radius_airgap_stator = radius_cylinder + airgap_stator_width;
// radius_airgap_rotor = radius_cylinder + airgap_stator_width + airgap_rotor_width;
// lc_cylinder = 0.1;
// lc_airgap = 0.05;

Include "parameters.geo";

Point(34) = {0, 0, 0, 1};
Point(35) = {0, radius_cylinder, 0, lc_cylinder};  
Point(36) = {radius_cylinder, 0, 0, lc_cylinder};  
Point(37) = {0, -radius_cylinder, 0, lc_cylinder}; 
Point(38) = {-radius_cylinder, 0, 0, lc_cylinder}; 

Point(26) = {0, radius_airgap_stator, 0, lc_airgap};
Point(27) = {radius_airgap_stator, 0, 0, lc_airgap};
Point(28) = {0, -radius_airgap_stator, 0, lc_airgap};
Point(29) = {-radius_airgap_stator, 0, 0, lc_airgap};//+
//+
Circle(1) = {38, 34, 35};
//+
Circle(2) = {35, 34, 36};
//+
Circle(3) = {36, 34, 37};
//+
Circle(4) = {37, 34, 38};
//+
Circle(5) = {29, 34, 26};
//+
Circle(6) = {26, 34, 27};
//+
Circle(7) = {27, 34, 28};
//+
Circle(8) = {28, 34, 29};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(6) = {2};

Mesh 2;
Physical Surface("air", 1) = {1};
//+
Physical Surface("cylinder", 2) = {6};
