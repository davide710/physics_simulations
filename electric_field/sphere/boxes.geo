// Gmsh project created on Fri Apr 26 11:27:31 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+

//+
Box(2) = {-5, -5, 0, 11, 11, -5};
//+
Box(3) = {-5, -5, 1, 11, 11, 5};
//+
Box(4) = {-5, -5, 0, 11, 5, 1};
//+
Box(5) = {-5, 1, 0, 11, 5, 1};
//+
Box(6) = {-5, 0, 0, 5, 1, 1};
//+
Box(7) = {1, 0, 0, 5, 1, 1};
//+
Physical Volume("box") = {1};
//+
Physical Volume("vacuum") = {4, 3, 6, 5, 2, 7};
