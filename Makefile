all: ga

ga:
	g++ -lm -o ga ./ga.o ./lk.o ./edgetour.o ./tsplib_io.o -std=c++11 -O3

run: ga
	./ga < ./input/cycle.in.$(n) 

clean:
	rm ga

ga.o: ga.cpp ga.h
	g++ -c ga.cpp -std=c++11 -lm -O3

lk.o: lk.cc lk.h
	g++ -c lk.cc -std=c++11 -lm -O3 -Wno-unused-result

edgetour.o: edgetour.cc tour.h
	g++ -c edgetour.cc -std=c++11 -lm -O3 -Wno-unused-result

tsplib_io.o: tsplib_io.cc tsplib_io.h
	g++ -c tsplib_io.cc -std=c++11 -lm -O3 -Wno-unused-result