Point(1) = {0, 0, 0};
//+
Point(2) = {10, 0, 0};
//+
Point(3) = {10, 10, 0};
//+
Point(4) = {0, 10, 0};
//+
Point(5) = {0, 0, 0.1};
//+
Point(6) = {10, 0, 0.1};
//+
Point(7) = {10, 10, 0.1};
//+
Point(8) = {0, 10, 0.1};
//+
Line(1) = {1, 2};
//+
Line(5) = {5, 6};
//+
Line(2) = {2, 3};
//+
Line(6) = {6, 7};
//+
Line(3) = {3, 4};
//+
Line(7) = {7, 8};
//+
Line(4) = {4, 1};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {2, 11, -6, -10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 10, -5, -9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -9, -4, 12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, -12, -3, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, 8, 5, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {3, 4, 1, 2};
//+
Plane Surface(6) = {6};
//+
Physical Surface("c1") = {5, 6, 1, 2, 3, 4};
//+
Surface Loop(1) = {5, 4, 3, 2, 6, 1};
//+
Volume(1) = {1};
//+
Physical Volume("c1_vol") = {1};
//+
Point(9) = {0, 0, 10.3};
//+
Point(10) = {10, 0, 10.3};
//+
Point(11) = {10, 10, 10.3};
//+
Point(12) = {0, 10, 10.3};
//+
Point(13) = {0, 0, 10.4};
//+
Point(14) = {10, 0, 10.4};
//+
Point(15) = {10, 10, 10.4};
//+
Point(16) = {0, 10, 10.4};
//+
Line(13) = {9, 10};
//+
Line(17) = {13, 14};
//+
Line(14) = {10, 11};
//+
Line(18) = {14, 15};
//+
Line(15) = {11, 12};
//+
Line(19) = {15, 16};
//+
Line(16) = {12, 9};
//+
Line(20) = {16, 13};
//+
Line(21) = {9, 13};
//+
Line(22) = {10, 14};
//+
Line(23) = {11, 15};
//+
Line(24) = {12, 16};
//+
//+
Curve Loop(7) = {20, 17, 18, 19};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {16, 13, 14, 15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {19, -24, -15, 23};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {13, 22, -17, -21};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {20, -21, -16, 24};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {18, -23, -14, 22};
//+
Plane Surface(12) = {12};
//+
Surface Loop(2) = {7, 11, 10, 8, 12, 9};
//+
Volume(2) = {2};
//+
Physical Surface("c2") = {7, 8, 11, 10, 12, 9};
//+
Physical Volume("c2_vol") = {2};
//+
//+
SetFactory("OpenCASCADE");
Sphere(3) = {5, 5, 10.4, 10, 0, Pi/2, 2*Pi};
//+
Sphere(4) = {5, 5, 0, 10, -Pi/2, 0, 2*Pi};
//+
Box(5) = {0, 0, 0.1, 10, 10, 10.2};
