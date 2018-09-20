FORT = ifort
##OPT = -convert big_endian -mcmodel=medium -shared-intel
##OPT = -mcmodel=large -shared-intel -heap-arrays
OPT = -heap-arrays
DBG = -traceback -check all

LIBS =  -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -lm
SRC = main.f90 bline.f90 sint.f90 blplot.f90 rotb_at_grid.f90
#SRC = main.f90 bline.f90 sint.f90
INC = paramete.f90 common.f90
OBJ = $(SRC:%.f90=%.o)
EXEC = bline_xy

main : $(OBJ) 
#	$(FORT) $(OPT) $(LIBS) $(OBJ) -o $(EXEC)
#	$(FORT) $(OPT) $(DBG) $(OBJ) -o $(EXEC)
	$(FORT) $(OPT) $(OBJ) -o $(EXEC)


%.o : %.f90 $(INC)
#	$(FORT) -c $< $(OPT) $(DBG)
	$(FORT) -c $< $(OPT)

clean:
	rm -f $(OBJ) $(EXEC)

