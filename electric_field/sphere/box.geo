//+
Point(3) = {-5, -5, -5, 1.0};
//+
Point(4) = {-5, -5, 5, 1.0};
//+
Point(5) = {-5, 5, -5, 1.0};
//+
Point(6) = {-5, 5, 5, 1.0};
//+
Point(7) = {5, -5, -5, 1.0};
//+
Point(8) = {5, -5, 5, 1.0};
//+
Point(9) = {5, 5, -5, 1.0};
//+
Point(10) = {5, 5, 5, 1.0};
//+
Line(4) = {5, 9};
//+
Line(5) = {9, 10};
//+
Line(6) = {10, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {3, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 4};
//+
Line(11) = {4, 3};
//+
Line(12) = {6, 4};
//+
Line(13) = {5, 3};
//+
Line(14) = {10, 8};
//+
Line(15) = {9, 7};
//+
Curve Loop(2) = {-7, -4, -5, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 10, 11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 15, -8, -13};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {6, 12, -10, -14};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, 13, -11, -12};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {5, 14, -9, -15};
//+
Plane Surface(7) = {7};
//+
Physical Surface("ext") = {2, 4, 7, 3, 6, 5};
//+
Surface Loop(2) = {2, 6, 4, 7, 5, 3};
//+
Surface Loop(3) = {1};
//+
Volume(2) = {2, 3};
//+
Physical Volume("mean") = {2, 3};