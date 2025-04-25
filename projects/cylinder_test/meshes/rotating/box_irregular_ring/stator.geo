Include "parameters.geo";

Point(34) = {0, 0, 0, 1};
Point(35) = {box_side / 2, box_side / 2, 0, lc_cube};  
Point(36) = {-box_side / 2, -box_side / 2, 0, lc_cube};  
Point(37) = {box_side / 2, -box_side / 2, 0, lc_cube}; 
Point(38) = {-box_side / 2, box_side / 2, 0, lc_cube}; 

Point(26) = {0, radius_airgap_stator, 0, lc_airgap};
Point(27) = {radius_airgap_stator, 0, 0, lc_airgap};
Point(28) = {0, -radius_airgap_stator, 0, lc_airgap};
Point(29) = {-radius_airgap_stator, 0, 0, lc_airgap};//+
//+
Line(1) = {38, 35};
//+
Line(2) = {35, 37};
//+
Line(3) = {37, 36};
//+
Line(4) = {36, 38};
//+
Circle(5) = {26, 34, 27};
//+
Circle(6) = {27, 34, 28};
//+
Circle(7) = {28, 34, 29};
//+
Circle(8) = {29, 34, 26};
//+
Curve Loop(1) = {8, 5, 6, 7};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(6) = {2};
//+
Physical Surface("air", 1) = {1};
//+
Physical Surface("cylinder", 2) = {6};

//Mesh 2; 
