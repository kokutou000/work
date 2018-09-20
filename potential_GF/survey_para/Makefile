FORT = ifort
FORTG = gfortran
DBGG = -fcheck=all -Wall
BIG = -convert big_endian

#LIBS =  -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -lm
DBG = -traceback -check all
OPT = -mcmodel=medium -shared-intel
OPTHEAP = -heap-arrays
SRC = parameter.f90 common.f90 valout.f90  main.f90 init.f90 vinit.f90 green_function.f90 diffs.f90
INC = parameter.f90 common.f90 valout.f90
OBJ = $(SRC:%.f90=%.o)
MOD = parameter.mod common.mod valout.mod
TAR = tmp
#EXEC = bline

### ifort
.SUFFIXES : .o .f90
.f90.o:
	${FORT} $(OPTHEAP) -c $<
${TAR} : ${OBJ}
	${FORT} $(OPTHEAP) -o $@ ${OBJ}

### gfortran 
#.SUFFIXES : .o .f90
#.f90.o:
#	${FORTG} $(DBGG) -c $<
#${TAR} : ${OBJ}
#	${FORTG} $(DBGG) -o $@ ${OBJ}
#
#main : $(OBJ) 
#	$(FORT) $(OPT) $(LIBS) $(OBJ)# -o $(EXEC)
#
#%.o : %.f90 $(INC)
#	$(FORT) -c $< $(OPT) 

clean:
	rm -f $(OBJ) *~

