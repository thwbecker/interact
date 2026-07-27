// splays example
// ==============
// gmsh -2 splays.geo
// gmsh -2 splays.geo -setnumber res_f 0.5
DefineConstant[ res = {100.0, Min 0, Max 10, Name "Domain resolution" } ];
DefineConstant[ res_f = {0.25, Min 0, Max 10, Name "Fault resolution" } ];
DefineConstant[ dip = {20, Min 0, Max 90, Name "Dipping angle" } ];
DefineConstant[ splay_dip1 = {50, Min 0, Max 90, Name "Splay dipping angle 1" } ];
DefineConstant[ splay_dip2 = {40, Min 0, Max 90, Name "Splay dipping angle 2" } ];
DefineConstant[ H = {60, Min 0, Max 10000, Name "Fault depth" } ];
DefineConstant[ S = {1000, Min 0, Max 10000, Name "Domain size" } ];

// some constants, dip in radiant
dip_rad = dip * Pi / 180.0;
splay_dip1_rad = splay_dip1 * Pi / 180.0;
splay_dip2_rad = splay_dip2 * Pi / 180.0;
dX = S;

// enable CAD modelling
SetFactory("OpenCASCADE");

X0 = -dX;
X1 = S * Cos(dip_rad) / Sin(dip_rad) + dX;
Y0 = -S;


Point(1) = {0.0, 0.0, 0.0, res_f};
Point(2) = {H * Cos(dip_rad) / Sin(dip_rad), -H, 0.0, res_f};
Point(3) = {S * Cos(dip_rad) / Sin(dip_rad), -S, 0.0, res};

// Create main fault
main_fault = newl; Line(main_fault) = {1, 2};
fault_extension = newl; Line(fault_extension) = {2, 3};

// Create splay faults with macro
Macro Splay
    det = Sin(dip_rad)*Cos(splay_dip2_rad) - Cos(dip_rad)*Sin(splay_dip2_rad);
    l2 = -Sin(dip_rad) * off / det;
    f1 = 0.3;
    f2 = 7.0;
    p1 = newp; Point(p1) = {off, 0.0, 0.0, res_f};
    p2 = newp; Point(p2) = {off + f1 * l2 * Cos(splay_dip1_rad),
                           -f1 * l2 * Sin(splay_dip1_rad), 0.0, res_f};
    p3 = newp; Point(p3) = {off + (l2 - f2) * Cos(splay_dip2_rad),
                           -(l2 - f2) * Sin(splay_dip2_rad), 0.0, res_f};
    Spline(splay_fault) = {p1, p2, p3};
Return

splay_fault1 = newl;
off = 30.0;
splay_fault = splay_fault1;
Call Splay;

splay_fault2 = newl;
off = 50.0;
splay_fault = splay_fault2;
Call Splay;

splay_fault3 = newl;
off = 70.0;
splay_fault = splay_fault3;
Call Splay;

splay_fault4 = newl;
off = 90.0;
splay_fault = splay_fault4;
Call Splay;

/*// Dirichlet "fault"*/
/*diri = news;*/
/*Line(diri) = {2, 3};*/

box = news;
Rectangle(box) = {X0, Y0, 0.0, X1-X0, -Y0};

// Intersect domain with fault
v() = BooleanFragments{ Surface{box}; Delete; }{ Line{main_fault, fault_extension, splay_fault1, splay_fault2, splay_fault3, splay_fault4}; Delete; };

// Recover lines with bounding boxes
eps = 1e-3;
top() = Curve In BoundingBox{X0-eps, -eps, -eps, X1+eps, eps, eps};
bottom() = Curve In BoundingBox{X0-eps, Y0-eps, -eps, X1+eps, Y0+eps, eps};
left() = Curve In BoundingBox{X0-eps, Y0-eps, -eps, X0+eps, eps, eps};
right() = Curve In BoundingBox{X1-eps, Y0-eps, -eps, X1+eps, eps, eps};
main_fault() = Curve In BoundingBox{-eps, Y0+eps, -eps, X1-eps, eps, eps};
main_fault() -= top();
main_fault() -= bottom();


// set mesh resolution
MeshSize{ PointsOf{Surface{:};} } = res;
MeshSize{ PointsOf{Line{main_fault(), splay_fault1, splay_fault2, splay_fault3};} } = res_f;

// 1 = free surface
Physical Curve(1) = {bottom(),top()};
// 3 = fault
Physical Curve(3) = {main_fault(), fault_extension, splay_fault1, splay_fault2, splay_fault3, splay_fault4};
// 5 = dirichlet
Physical Curve(5) = {left(),right()};
Physical Surface(1) = {v()};

Mesh.MshFileVersion = 2.2;
