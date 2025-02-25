//+
SetFactory("OpenCASCADE");
h1 = 1;
h2 = 1;
h = 0.32;
l = 0.32;
r = 0.02;

Point(1) = {0, 0, 0, h1};
Point(2) = {l, 0, 0, h1};
Point(3) = {l, h, 0, h1};
Point(4) = {0, h, 0, h1};

//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(6) = {0.16, 0.16, 0, r, 0, 2*Pi};  // Center of the circle and radius

//+
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {6};
Plane Surface(1) = {1, 2};

