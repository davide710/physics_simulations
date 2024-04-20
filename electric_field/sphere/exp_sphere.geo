// Gmsh project created on Sat Apr 20 15:26:21 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.1, 0, 2*Pi};
//+
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Extrude {{0, 0, 0}, {0, 0, 0}, Pi} {
  Point{1}; Curve{1}; Surface{1}; 
}
