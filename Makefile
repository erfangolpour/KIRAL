all: debug

debug:
	g++ -g -Wall -o main main.cpp -I. -L./minimap2 -lminimap2 -lz

release:
	g++ -O3 -o main main.cpp -I. -L./minimap2 -lminimap2 -lz

clean:
	rm -f main main.o