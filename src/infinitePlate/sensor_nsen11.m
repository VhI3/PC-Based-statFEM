% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)

%
%  Matlab mesh
% square_plate_with_hole2, Created by Gmsh 3.0
% ASCII
clear sensor;
sensor.nbNod = 11;
sensor.POS = [
0 0 0;
0.32 0 0;
0.32 0.32 0;
0 0.32 0;
0.18 0.16 0;
0.1724697960371747 0.1756366296493606 0;
0.1555495813208737 0.1794985582436365 0;
0.1419806226419516 0.1686776747823512 0;
0.1419806226419516 0.1513223252176489 0;
0.1555495813208737 0.1405014417563635 0;
0.1724697960371747 0.1443633703506394 0;
];
sensor.MAX = max(sensor.POS);
sensor.MIN = min(sensor.POS);
sensor.LINES =[
 1 2 0
 2 3 0
 3 4 0
 4 1 0
 5 6 0
 6 7 0
 7 8 0
 8 9 0
 9 10 0
 10 11 0
 11 5 0
];
sensor.TRIANGLES =[
 9 1 10 0
 11 2 5 0
 7 4 8 0
];
sensor.QUADS =[
 4 1 9 8 0
 10 1 2 11 0
 3 4 7 6 0
 5 2 3 6 0
];
sensor.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
];
