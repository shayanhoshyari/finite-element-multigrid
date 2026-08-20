Point(1) = {0., 0., 0., 0.4};
Point(2) = {1., 0., 0., 0.4};
Point(3) = {1., 1., 0., 0.4};
Point(4) = {0., 1., 0., 0.4};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = { 1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(101) = { 1 };
Physical Line(101) = { 1 };
Physical Line(102) = { 2 };
Physical Line(103) = { 3 };
Physical Line(104) = { 4 };
