## Discrete Element Method (DEM) Simulations
The code implements a DEM simluation of grains in bi-periodic shear. Grain interactions are easily customisable in contact.cpp, where they currently include elasticity, viscous forces, cohesion and friction. 

The source c++ code was written by Pierre Rognon and modified by Prashida Kharel and Matthew Macaulay.
The post processing in Matlab was written by Matthew Macaulay.

## Installation and running
Open DEM_V2.pro using Qt Creator and press run.

Alternatively, once Qt Creator has built it, you can call the binary file directly from the command line interface. The binary is located at the <build dir>/DEM_V2.app/Contents/MacOS/DEM_V2

Once running, the program will ask for inputs using a keyword. For example using the action keyword: ACTION CREATE_RANDOM. 

# Computational Time
Lower shear rates take longer, above 10^-3 is reasonable.
The pairwise-distance computation is O(n^2). n=10^4 is reasonable.   


## Post Processing
In Matlab, create an instance of a e.g. PostP(dir) object with the same directory as the data <pp = PostP(/path/to/data);>. Then use any of the PostP functions on this e.g.  pp.makeMovie();