// Gmsh project created on Thu Apr 25 17:09:23 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 5, 0, 2*Pi};
//+
Extrude {0, 0, 10} {
  Point{2}; Curve{2}; Point{1}; Curve{1}; 
}
//+
Curve Loop(3) = {6};
//+
Curve Loop(4) = {8};
//+
Plane Surface(3) = {3, 4};
//+
Curve Loop(5) = {2};
//+
Curve Loop(6) = {1};
//+
Plane Surface(4) = {5, 6};
//+
Surface Loop(1) = {2, 3, 1, 4};
//+
Volume(1) = {1};
//+
Curve Loop(9) = {8};
//+
Plane Surface(7) = {9};
//+
Curve Loop(10) = {1};
//+
Plane Surface(8) = {10};
//+
Surface Loop(2) = {2, 7, 8};
//+
Volume(2) = {2};
//+
Physical Volume("vacuum") = {1};
//+
Physical Volume("wire") = {2};
