all: area area-mp

area: area2.c
	gcc -O3 area2.c -o area2

area-mp: area2.c
	gcc -O3 -fopenmp area2.c -o area2-mp
