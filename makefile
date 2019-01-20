OPTIONS = -w #-O0
COMP = gfortran

OUT = exe
PRGM = functions.f90 main.f90

exe : $(PRGM)
	$(COMP) $(OPTIONS) -o $(OUT) $(PRGM)

clean :
	rm *.o
	rm *.mod

plot :
	gnuplot -p plot.gnu


anim :
	gnuplot -p anim.gnu
