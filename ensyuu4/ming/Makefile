all: ming ming-mp ming-mp-bh

ming: ming.c mypng.o genpng.o judge.o
	gcc -O3 ming.c mypng.o genpng.o judge.o -o ming -lm -lpng

ming-mp: ming.c mypng.o genpng.o judge.o
	gcc -O3 -fopenmp ming.c mypng.o genpng.o judge.o -o ming-mp -lm -lpng

ming-mp-bh: ming-bh.c mypng.o genpng.o judge.o
	gcc -O3 -fopenmp ming-bh.c mypng.o genpng.o judge.o -o ming-mp-bh -lm -lpng

mypng.o: mypng.c mypng.h genpng.h xydata.h
	gcc -O3 mypng.c -c

genpng.o: genpng.c genpng.h
	gcc -O2 genpng.c -c

m1.png: ming Makefile
	./ming 20 1.07 0.7 0 > m1.png  # 2pi/6

m2.png: ming Makefile
	./ming 20 2.6 2.8 0  > m2.png  # pi - 2pi/12

otherPNGs: ming Makefile
	./ming 20 3.17 4.3 0 > a2.png  # 2pi/2
	./ming 20 2.14 4.3 0 > a3.png  # 2pi/3
	./ming 20 1.6 0.3 0 > a4.png   # 2pi/4
	./ming 20 1.27 0.3 0 > a5.png  # 2pi/5
	./ming 20 1.07 0.3 0 > a6.png  # 2pi/6
	./ming 20 0.8 0.3 0 > a8.png   # 2pi/8
	./ming 20 0.713 0.3 0 > a9.png # 2pi/9
	./ming 20 0.54 0.3 0 > a12.png # 2pi/12
