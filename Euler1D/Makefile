
all: gridFunc1D.o simulation.o Euler1D.o parameterReader.o mvector.o timeIntegrator.o
	g++ gridFunc1D.o simulation.o Euler1D.o parameterReader.o mvector.o timeIntegrator.o -o euler1D
	mv euler1D exe

gridFunc1D.o: gridFunc1D.cpp gridFunc1D.hpp
	g++ -g -c gridFunc1D.cpp -o gridFunc1D.o

Euler1D.o: Euler1D.cpp Euler1D.hpp
	g++ -g -c Euler1D.cpp -o Euler1D.o

simulation.o: simulation.cpp 
	g++ -g -c simulation.cpp -o simulation.o

parameterReader.o: parameterReader.cpp parameterReader.hpp
	g++ -g -c parameterReader.cpp -o parameterReader.o

mvector.o: mvector.hpp mvector.cpp
	g++ -g -c mvector.cpp -o mvector.o

timeIntegrator.o: timeIntegrator.cpp timeIntegrator.hpp
	g++ -g -c timeIntegrator.cpp -o timeIntegrator.o


clean:
	rm *.o
