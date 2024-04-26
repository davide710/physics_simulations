// Gmsh project created on Fri Apr 26 11:27:31 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Box(2) = {-5, -5, -5, 11, 11, 11};
//+
Surface Loop(3) = {12, 7, 9, 11, 10, 8};
//+
Surface Loop(4) = {6, 1, 3, 5, 4, 2};
//+
Volume(3) = {3, 4};
//+
Physical Volume("box") = {1};
//+
Physical Volume("vacuum") = {3};
