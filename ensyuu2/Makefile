all: area area-mp

area: area.c
	gcc -O3 area.c -o area

area-mp: area.c
	gcc -O3 -fopenmp area.c -o area-mp
