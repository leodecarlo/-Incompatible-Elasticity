//+
SetFactory("OpenCASCADE");
//Rectangle(1) = {0, 0, 0, 2.5, 0.41, 0};
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {1.0, 0.0, 0.0, 1.0};
Point(3) = {1.0, 1.0, 0.0, 1.0};
Point(4) = {0.0, 1.0, 0.0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};

//LH edge of circle
//Point(5) = {0.15, 0.2, 0, 0.1};
//C
//Point(6) = {0.2, 0.2, 0, 0.1};
//RH edge of circle
//Point(7) = {0.25, 0.2, 0, 0.1};
//Upper meeting pt of flag with pole: REFINE near here
//Point(8) = {0.248989795, 0.21, 0, 0.1};
//Lower meeting pt of flag with pole: REFINE near here
//Point(9) = {0.248989795, 0.19, 0, 0.1};
//Corner of the flag above A : REFINE near here
//Point(10) = {0.6, 0.21, 0, 0.1};
//Corner of the flag below A: REFINE near here
//Point(11) = {0.6, 0.19, 0, 0.1};
//Circle(5) = {5, 6, 7};
//Circle(6) = {7, 6, 5};
//Line(7) = {8, 10};
//Line(8) = {10, 11};
//Line(9) = {11, 9};
//Line Loop(2) = {5, 6};
//Line Loop(2) = {5, 6};
//Line Loop(3) = {7, 8, 9};
//Have to `delete' the point C?
//Delete{Point{6};}

Circle(5) = {0.5, 0.5, 0, 0.1, 0, 2*Pi};
Line Loop(2) = {5};

Plane Surface(1) = {1, 2};

// --- Define physical groups for boundaries ---
Physical Curve("Left")   = {4};
Physical Curve("Right")  = {2};
Physical Curve("Top")    = {3};
Physical Curve("Bottom") = {1};
Physical Curve("Hole")   = {5};

// --- Define physical surface for the plate ---
Physical Surface("Plate") = {1};


//gmsh -2 -clscale 0.1 -o perforated_plate_clscale0.1.msh perforated_plate.geo


