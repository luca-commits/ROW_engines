// scale = 1; 
//    x0_ = -25; y0_ = -20; w_ =- 20*scale ; h_ =  40*scale; r_ = -1.5*scale;

    
//     // Points for bottom side
//     p1 = newp; Point(p1) = {x0_ + r_, y0_ +r_, 0, 1.0};
//     p2 = newp; Point(p2) = {x0_ + w_ - r_, y0_ + r_, 0, 1.0};
    
//     // Bottom-right arc
//     pc1 = newp; Point(pc1) = {x0_ + w_ - r_, y0_ + r_, 0, 1.0};
//     p3 = newp; Point(p3) = {x0_ + w_, y0_ + r_, 0, 1.0};
    
//     // Right side
//     p4 = newp; Point(p4) = {x0_ + w_, y0_ + h_ - r_, 0, 1.0};
    
//     // Top-right arc
//     pc2 = newp; Point(pc2) = {x0_ + w_ - r_, y0_ + h_ - r_, 0, 1.0};
//     p5 = newp; Point(p5) = {x0_ + w_ - r_, y0_ + h_, 0, 1.0};
    
//     // Top side
//     p6 = newp; Point(p6) = {x0_ + r_, y0_ + h_, 0, 1.0};
    
//     // Top-left arc
//     pc3 = newp; Point(pc3) = {x0_ + r_, y0_ + h_ - r_, 0, 1.0};
//     p7 = newp; Point(p7) = {x0_, y0_ + h_ - r_, 0, 1.0};
    
//     // Left side
//     p8 = newp; Point(p8) = {x0_, y0_ - r_, 0, 1.0};
    
//     // Bottom-left arc
//     pc4 = newp; Point(pc4) = {x0_ + r_, y0_ + r_, 0, 1.0};
//     p9 = newp; Point(p9) = {x0_ + r_, y0_, 0, 1.0};
    
//     // Create curves
//     l1 = newl; Line(l1) = {p1, p2};
//     a1 = newl; Circle(a1) = {p2, pc1, p3};
//     l2 = newl; Line(l2) = {p3, p4};
//     a2 = newl; Circle(a2) = {p4, pc2, p5};
//     l3 = newl; Line(l3) = {p5, p6};
//     a3 = newl; Circle(a3) = {p6, pc3, p7};
//     l4 = newl; Line(l4) = {p7, p8};
//     a4 = newl; Circle(a4) = {p8, pc4, p1};
    
//     // Create curve loop and surface
//     cl = newl; Curve Loop(cl) = {l1, a1, l2, a2, l3, a3, l4, a4};
//     surf = news; Plane Surface(surf) = {cl};


p1 = newp; Point(p1 = )
