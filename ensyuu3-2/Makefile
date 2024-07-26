all: area2 area2-mp sum sum-mp 

area2: area2.c
	gcc -O3 area2.c -o area2

area2-mp: area2.c
	gcc -O3 -fopenmp area2.c -o area2-mp

sum: sum.c
	gcc -O3 sum.c -o sum

sum-mp: sum.c
	gcc -O3 -fopenmp sum.c -o sum-mp

