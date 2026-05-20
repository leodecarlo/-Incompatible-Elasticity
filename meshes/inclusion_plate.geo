//+
SetFactory("OpenCASCADE");
//Rectangle(1) = {0, 0, 0, 2.5, 0.41, 0};
coarse_length = 1.0;
Point(1) = {0.0, 0.0, 0.0, coarse_length};
Point(2) = {1.0, 0.0, 0.0, coarse_length};
Point(3) = {1.0, 1.0, 0.0, coarse_length};
Point(4) = {0.0, 1.0, 0.0, coarse_length};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Line Loop(1) = {1, 2, 3, 4};
Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

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

fine_length = 1.0;
//Circle(5) = {0.5, 0.5, 0, 0.1, 0, 2*Pi};
Point(5) = {0.5, 0.5, 0, fine_length};
Point(6) = {0.5, 0.65, 0, fine_length};
Point(7) = {0.65, 0.5, 0, fine_length};
Point(8) = {0.5, 0.35, 0, fine_length};
Point(9) = {0.35, 0.5, 0, fine_length};
Circle(5) = {6, 5, 9};
Circle(6) = {9, 5, 8};
Circle(7) = {8, 5, 7};
Circle(8) = {6, 5, 7};
//Line Loop(2) = {5};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(2) = {2};
//Plane Surface(2) = {5}; // assuming I can use the Circle as a Curve Loop?

// Define the outer region by "cutting" the circle out of the square
// The outer surface is a combination of the square's loop and the circle's loop
Plane Surface(3) = {1, 2};

// Define physical groups for the different regions
Physical Surface("Circle") = {2};
Physical Surface("Square") = {3};


// From Google AI:
//// Create a mesh size field for the subregion
//// We use a Distance field to control mesh size based on proximity to the refined surface.
//// This ensures a smooth transition between the coarse and fine meshes.
//Field[1] = MathEval;
//Field[1].F = "F2/50.0 + (1-F2)*(1-F2)*coarse_length";
//// Create a second MathEval field for the finer region to be more explicit.
//Field[2] = MathEval;
//Field[2].F = "fine_length";
//// Combine the fields
//Field[3] = Box;
//Field[3].VIn = 2; // Use Field[2] inside the box
//Field[3].VOut = 1; // Use Field[1] outside the box
////Field[3].XMin = 4; // FA: not sure what these lines were for
////Field[3].XMax = 6;
////Field[3].YMin = 4;
////Field[3].YMax = 6;
//// Specify the final mesh size field to be the minimum of all defined fields
//Field[4] = Min;
//Field[4].FieldsList = {1, 3};
//// Use the combined field as the background mesh size
//Background Field = 4;
//// Generate the 2D mesh
////Mesh.Algorithm = 6; // Choose a 2D unstructured algorithm (e.g., MeshAdapt)
////Mesh 2;

//gmsh -2 -clscale 0.04 -o inclusion_plate_clscale0.04.msh inclusion_plate.geo


