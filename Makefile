all: ga

ga: ga.cpp
	g++ -o ga -O3 ga.cpp

run: ga
	./ga < ./cycle.in > ./cycle.out

clean:
	rm ga
