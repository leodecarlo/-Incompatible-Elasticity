//----------------------------------------------------
// Perforated Plate Geometry for Firedrake Simulation
// Square domain with a circular hole in the center
// Compatible with Gmsh 4.x and Firedrake
//----------------------------------------------------

SetFactory("OpenCASCADE");

// === Outer square points ===
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {1.0, 0.0, 0.0, 1.0};
Point(3) = {1.0, 1.0, 0.0, 1.0};
Point(4) = {0.0, 1.0, 0.0, 1.0};

// === Outer boundary lines ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};

// === Inner circular hole ===
// Define center and four quarter points
Point(5) = {0.5, 0.6, 0, 1.0};
Point(6) = {0.6, 0.5, 0, 1.0};
Point(7) = {0.5, 0.4, 0, 1.0};
Point(8) = {0.4, 0.5, 0, 1.0};
Point(9) = {0.5, 0.5, 0, 1.0}; // center of circle

// Define four arcs around the circle
Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};

// Combine arcs into a circular loop
Line Loop(2) = {5, 6, 7, 8};

// === Define surface (outer minus inner) ===
Plane Surface(1) = {1, 2};

// === Define physical boundaries ===
// Label each edge and the hole boundary for Firedrake BCs
Physical Curve("Left")   = {4};
Physical Curve("Right")  = {2};
Physical Curve("Top")    = {3};
Physical Curve("Bottom") = {1};
Physical Curve("Hole")   = {5, 6, 7, 8};

// Define physical surface (the solid plate)
Physical Surface("Plate") = {1};

// === Notes ===
// Use command below to generate the mesh file
 gmsh -2 -clscale 0.1 -o perforated_plate_clscale0.1.msh perforated_plate.geo
//----------------------------------------------------

//gmsh -2 -clscale 0.1 -o perforated_plate_clscale0.1.msh perforated_plate.geo


