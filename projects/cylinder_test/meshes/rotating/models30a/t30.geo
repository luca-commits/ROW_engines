Include "t30_data.geo" ;

Solver.AutoShowLastStep = 1;
Mesh.Algorithm = 1 ;

// -------------------------------------------------------
// Some characteristic lengths
//--------------------------------------------------------

lr1 = r1/15 ;
lr2 = r1/45 ;
lr3 = lr2 ;
lr4 = r4/10 ;
lr5 = lr4 ;

lc = lr1 ;

//--------------------------------------------------------

cen = 1;
Point(cen) = {0,0,0,lc};

Printf("ROTOR_FE: %d A",ROTOR_FE);

Include "t30_rotor.geo";
Include "t30_stator.geo";


// Hide { Point{ Point '*' }; }
// //Hide { Line{ Line '*' }; }
// Show { Line{ nicepos_rotor[], nicepos_stator[] }; }

// //For post-processing...
// //View[0].Light = 0;
// View[0].NbIso = 25; // Number of intervals
// View[0].IntervalsType = 1;
// //+
// Show "*";
// //+
// Hide {
//   Point{1}; Point{2}; Point{3}; Point{4}; Point{5}; Point{6}; Point{7}; Point{8}; Point{9}; Point{10}; Point{11}; Point{12}; Point{13}; Point{14}; Point{15}; Point{16}; Point{17}; Point{18}; Point{19}; Point{22}; Point{25}; Point{26}; Point{27}; Point{28}; Point{29}; Point{30}; Point{31}; Point{32}; Point{33}; 
// }

// //+
// Plane Surface(1040) = {15, 1038};
// Mesh 2; 