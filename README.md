# The Kepler DR24 Robovetter

This code will run the Kepler Q1-Q17 DR24 Robovetter on the inputs given in the Q1-Q17 DR24 KOI catalog paper (Coughlin et al. 2015) and produce outputs given in the catalog paper and at the NASA Exoplanet Q1-Q17 DR24 KOI table. This will work for both the real dataset, as well as the dataset that contains artificially injected trasniting planets.

To compile the code, either type "make" to use the Makefile, or compile via:

g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp

To run the code, supply it an input file and an output filename. For example:

./robovet RoboVetter-Input.txt RoboVetter-Output.txt

for the real data, or
 
./robovet RoboVetter-Inject-Input.txt RoboVetter-Inject-Output.txt

for the artifically injected transit data.

Note that in this case, "RoboVetter-Input.txt" corresponds to Table 3 of the paper, "RoboVetter-Output.txt" corresponds to Table 5, and "RoboVetter-Inject-Input.txt" and "RoboVetter-Inject-Output.txt" corresponds to Table 6.
