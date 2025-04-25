Include "parameters.geo";


Point(1) = {0, 0, 0, 1};
Point(2) = {0, radius_ring_internal, 0, lc_ring_internal};  
Point(3) = {radius_ring_internal, 0, 0, lc_ring_internal};  
Point(4) = {0, -radius_ring_internal, 0, lc_ring_internal}; 
Point(5) = {-radius_ring_internal, 0, 0, lc_ring_internal}; 

Point(6) = {0, radius_ring_external, 0, lc_ring_external};  
Point(7) = {radius_ring_external, 0, 0, lc_ring_external};  
Point(8) = {0, -radius_ring_external, 0, lc_ring_external}; 
Point(9) = {-radius_ring_external, 0, 0, lc_ring_external}; 

// //right irregularity
box_width = 0.1;
box_heidth = 0.2;
lc_box = 0.1;

Point(10) = {radius_ring_external, -box_heidth/2, 0, lc_box};
Point(11) = {radius_ring_external + box_width, -box_heidth/2, 0, lc_box};
Point(12) = {radius_ring_external + box_width, +box_heidth/2, 0, lc_box};
Point(13) = {radius_ring_external, box_heidth/2, 0, lc_box};

Point(14) = {radius_ring_internal, -box_heidth/2, 0, lc_box};
Point(15) = {radius_ring_internal, box_heidth/2, 0, lc_box};
Point(16) = {radius_ring_internal - box_width, +box_heidth/2, 0, lc_box};
Point(17) = {radius_ring_internal - box_width, -box_heidth/2, 0, lc_box};

Point(18) = {-radius_ring_internal, -box_heidth/2, 0, lc_box};
Point(19) = {-radius_ring_internal + box_width, -box_heidth/2, 0, lc_box};
Point(20) = {-radius_ring_internal + box_width, +box_heidth/2, 0, lc_box};
Point(21) = {-radius_ring_internal, box_heidth/2, 0, lc_box};

Point(22) = {-radius_ring_external, -box_heidth/2, 0, lc_box};
Point(23) = {-radius_ring_external, box_heidth/2, 0, lc_box};
Point(24) = {-radius_ring_external - box_width, +box_heidth/2, 0, lc_box};
Point(25) = {-radius_ring_external - box_width, -box_heidth/2, 0, lc_box};

lc_airbox = 0.1;
radius_airbox = 3;
Point(30) = {0, radius_airbox, 0, lc_airbox};  
Point(31) = {radius_airbox, 0, 0, lc_airbox};  
Point(32) = {0, -radius_airbox, 0, lc_airbox}; 
Point(33) = {-radius_airbox, 0, 0, lc_airbox}; 


// Define points and curves for outer circle
Point(39) = {0, radius_airgap_rotor, 0, lc_airgap};
Point(40) = {radius_airgap_rotor, 0, 0, lc_airgap};
Point(41) = {0, -radius_airgap_rotor, 0, lc_airgap};
Point(42) = {-radius_airgap_rotor, 0, 0, lc_airgap};




//+
Circle(1) = {39, 1, 40};
//+
Circle(2) = {40, 1, 41};

//+
Circle(3) = {41, 1, 42};
//+
Circle(4) = {42, 1, 39};
//+
Circle(5) = {14, 1, 4};
//+
Circle(6) = {4, 1, 18};
//+
Circle(7) = {21, 1, 2};
//+
Circle(8) = {2, 1, 15};
//+
Circle(9) = {6, 1, 13};
//+
Circle(10) = {10, 1, 8};
//+
Circle(11) = {8, 1, 22};
//+
Circle(12) = {23, 1, 6};
//+
Circle(13) = {30, 1, 31};
//+
Circle(14) = {31, 1, 32};
//+
Circle(15) = {32, 1, 33};
//+
Circle(16) = {33, 1, 30};
//+
Line(17) = {21, 18};
//+
Line(18) = {21, 20};
//+
Line(19) = {20, 19};
//+
Line(20) = {19, 18};
//+
Line(21) = {22, 25};
//+
Line(22) = {24, 24};
//+
Line(23) = {24, 25};
//+
Line(24) = {23, 24};
//+
Line(25) = {23, 22};
//+
Line(26) = {16, 17};
//+
Line(27) = {17, 14};
//+
Line(28) = {15, 14};
//+
Line(29) = {15, 16};
//+
Line(30) = {10, 11};
//+
Line(31) = {12, 11};
//+
Line(32) = {12, 13};
//+
Line(33) = {13, 10};
//+
Curve Loop(1) = {16, 13, 14, 15};
//+
Curve Loop(2) = {12, 9, -32, 31, -30, 10, 11, 21, -23, -24};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {9, 33, 10, 11, -25, 12};
//+
Curve Loop(4) = {17, -6, -5, -28, -8, -7};
//+
Plane Surface(2) = {3, 4};
//+
Curve Loop(5) = {7, 8, 29, 26, 27, 5, 6, -20, -19, -18};
//+
Curve Loop(6) = {2, 3, 4, 1};
//+
Plane Surface(3) = {5, 6};

Physical Surface("air", 1) = {1, 3};
//+
Curve Loop(7) = {23, -21, -25, 24};
//+
Plane Surface(4) = {7};
//+
Curve Loop(8) = {18, 19, 20, -17};
//+
Plane Surface(5) = {8};
//+
Curve Loop(9) = {29, 26, 27, -28};
//+
Plane Surface(8) = {9};
//+
Curve Loop(10) = {32, 33, 30, -31};
//+
Plane Surface(7) = {10};
//+
Physical Surface("ring", 3) = {2, 4, 5, 8, 7};
Mesh 2; //+//+
