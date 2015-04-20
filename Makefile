all: ga

ga: ga.cpp
	g++ -o ga -O3 ga.cpp

run: ga
	./ga < ./input/cycle.in.$(n) 

clean:
	rm ga
