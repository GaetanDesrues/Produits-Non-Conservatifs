OPTIONS = -w #-O0
COMP = gfortran

OUT = exe
PRGM = test.f90 main.f90
# PRGM = aulerMieussens.f90

exe : $(PRGM)
	$(COMP) $(OPTIONS) -o $(OUT) $(PRGM)

clean :
	rm *.o
	rm *.mod

cleanoutput :
	rm Output/*.txt

plot :
	gnuplot -p plot.gnu

anim :
	gnuplot -p anim.gnu
