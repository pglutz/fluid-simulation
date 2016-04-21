vis: Particles.h
	g++ -std=c++0x main.cpp Particles.cpp -o sim -lGL -lGLU -lglut -lopencv_core -lopencv_highgui -I.
