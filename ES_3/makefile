main : main.cpp error.o random.o
	g++ -o main main.cpp error.o random.o
RNGExtend.o : RNGExtend.cpp RNGExtend.h
	g++ -c -o RNGExtend.o RNGExtend.cpp 
error.o : error.cpp error.h
	g++ -c -o error.o error.cpp 
clear : 
	rm *.o