SetFactory("OpenCASCADE");
Merge "simulation1.step";
//+
Physical Surface("ice", 19) = {2};
//+
Physical Surface("saturated_shale_and_silt", 20) = {1, 3};
//+
Physical Surface("gneiss", 21) = {4};
//+
Physical Curve("boundary", 22) = {6, 5, 4, 3, 13, 10, 9, 12};
