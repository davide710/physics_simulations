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
Surface Loop(1) = {5, 4, 3, 2, 6, 1};
//+
Volume(1) = {1};
//+
Point(9) = {0, 0, 0.3};
//+
Point(10) = {10, 0, 0.3};
//+
Point(11) = {10, 10, 0.3};
//+
Point(12) = {0, 10, 0.3};
//+
Point(13) = {0, 0, 0.4};
//+
Point(14) = {10, 0, 0.4};
//+
Point(15) = {10, 10, 0.4};
//+
Point(16) = {0, 10, 0.4};
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
Physical Surface("c+_sup") = {8};
//+
Physical Surface("c-_sup") = {5};
//+
Physical Volume("c+_vol") = {2};
//+
Physical Volume("c-_vol") = {1};
//+
SetFactory("OpenCASCADE");
Box(3) = {0, 0, 0.1, 10, 10, 0.2};
//+
Box(4) = {0, 0, 0, -5, 15, 0.4};
//+
Box(5) = {0, 10, 0, 15, 5, 0.4};
//+
Box(6) = {10, 10, 0, 5, -15, 0.4};
//+
Box(7) = {10, 0, 0, -15, -5, 0.4};
//+
Box(8) = {-5, -5, 0, 20, 20, -5};
//+
Box(9) = {-5, -5, 0.4, 20, 20, 5};
//+
Physical Surface("ext") = {54, 49, 52, 50, 51, 43, 46, 44, 45, 47, 39, 33, 37, 19, 22, 28, 26, 32};
//+
Physical Volume("vacuum") = {9, 8, 3, 4, 7, 6, 5};
