all: area area-mp area2 area2-mp

area: area.c
	gcc -O3 area.c -o area

area-mp: area.c
	gcc -O3 -fopenmp area.c -o area-mp

area2: area2.c
	gcc -O3 area2.c -o area2

area2-mp: area2.c
	gcc -O3 -fopenmp area2.c -o area2-mp
