# The Kepler DR24 Robovetter

This code will run the Kepler Q1-Q17 DR24 Robovetter on the inputs given in the Q1-Q17 DR24 KOI catalog paper (Coughlin et al. 2015) and produce outputs given in the catalog paper and at the NASA Exoplanet Q1-Q17 DR24 KOI table. This will work for both the real dataset, as well as the dataset that contains artificially injected trasniting planets.

To compile the code, either type "make" to use the Makefile, or compile via:

g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp

To run the code, supply it an input file and an output filename. For example:

g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp

for the real data, or:

g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp

for the artificially injected transits.


