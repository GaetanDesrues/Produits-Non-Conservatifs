OPTIONS = -w #-O0
COMP = gfortran

OUT = exe

exe : test.f90 main.f90
	$(COMP) $(OPTIONS) -o $(OUT) test.f90 main.f90

clean :
	rm *.o
	rm *.mod

cleanoutput :
	rm Output/*.txt

plot :
	gnuplot -p plot.gnu

anim :
	gnuplot -p anim.gnu
