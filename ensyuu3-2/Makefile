all: area area-mp area3s area3s-mp area3r area3r-mp

area: area3d.c
	gcc -O3 area3d.c -o area3d

area-mp: area3d.c
	gcc -O3 -fopenmp area3d.c -o area3d-mp

area3s:	area3s.c
	gcc -O3 area3s.c -o area3s

area3s-mp:	area3s.c
	gcc -O3 -fopenmp area3s.c -o area3s-mp

area3r:	area3r.c
	gcc -O3 area3r.c -o area3r

area3r-mp: area3r.c
	gcc -O3 -fopenmp area3r.c -o area3r-mp

	
