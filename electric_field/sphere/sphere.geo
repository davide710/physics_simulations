// Gmsh project created on Fri Apr 19 18:24:20 2024
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0.000, 0.000, 0.000, 0.1, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("sphere_sup") = {1};
//+
Physical Volume("sphere") = {1};
//+
Sphere(2) = {0.000, 0.000, 0.000, 5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("large_sphere_sup") = {2};
//+
Physical Volume("large_sphere") = {2, 1};
//+
