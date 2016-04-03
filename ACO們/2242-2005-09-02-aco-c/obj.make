CC = g++
objfn : test.o  obj.o
	$(CC)  -g -o  objfn test.o obj.o -lm 
test.o : test.cpp
	$(CC)  -g -c test.cpp
obj.o:obj.c parameters.h 
	g++  -g -c obj.c 
	
