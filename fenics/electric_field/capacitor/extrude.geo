// Gmsh project created on Wed Apr 24 17:16:15 2024
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 11, 0, 1.0};
//+
Point(4) = {0, 11, 0, 1.0};
//+
Point(5) = {5, 5, 0, 1.0};
//+
Point(6) = {15, 5, 0, 1.0};
//+
Point(7) = {15, 5.2, 0, 1.0};
//+
Point(8) = {5, 5.2, 0, 1.0};
//+
Point(9) = {5, 5.8, 0, 1.0};
//+
Point(10) = {15, 5.8, 0, 1.0};
//+
Point(11) = {15, 6, 0, 1.0};
//+
Point(12) = {5, 6, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 9};
//+
Curve Loop(1) = {11, 12, 9, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, 1, 2};
//+
Plane Surface(3) = {1, 2, 3};
//+
Physical Surface("c+") = {1};
//+
Physical Surface("c-") = {2};
//+
Physical Surface("vacuum") = {3};
//+
Extrude {0, 0, 5} {
  Point{12}; Point{11}; Point{10}; Point{9}; Point{8}; Point{7}; Point{6}; Point{5}; Point{1}; Point{4}; Point{3}; Point{2}; Curve{2}; Curve{3}; Curve{4}; Curve{1}; Curve{5}; Curve{6}; Curve{7}; Curve{8}; Curve{9}; Curve{10}; Curve{11}; Curve{12}; Surface{2}; Surface{1}; Surface{3}; 
}
//+
Physical Volume("c+_vol") = {2};
//+
Physical Volume("c-_vol") = {1};
//+
Physical Volume("vacuum_vol") = {3};
