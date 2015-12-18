all : robovetter

robovetter : DR24-RoboVetter.cpp
	g++ -std=c++11 -O3 -o robovet DR24-RoboVetter.cpp
