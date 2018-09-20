!FC = f90
FC = sxmpif90

EXEC = nl3d_mpi02k_032x032x0512
OBJ = constants.o common.o mpiut.o job.o pset.o iset.o rkg.o out.o mhd.o eflux.o urand1.o
INC = constants.f90 

OPT = -ftrace
#OPT = -P auto -ftrace

all : main.f90 $(OBJ)
	$(FC) $(OPT) -o $(EXEC) main.f90 $(OBJ)

constants.o : constants.f90
	$(FC) $(OPT) -c constants.f90

common.o : constants.o common.f90
	$(FC) $(OPT) -c common.f90

mpiut.o : common.o mpiut.f90
	$(FC) $(OPT) -c mpiut.f90

pset.o : common.o mpiut.o pset.f90
	$(FC) $(OPT) -c pset.f90

job.o : constants.o common.f90 job.f90
	$(FC) $(OPT) -c job.f90

iset.o : common.o iset.f90
	$(FC) $(OPT) -c iset.f90

iset_model.o : common.o iset_model.f90
	$(FC) $(OPT) -c iset_model.f90

rkg.o : common.o rkg.f90
	$(FC) $(OPT) -c rkg.f90

out.o : common.o out.f90
	$(FC) $(OPT) -c out.f90

eflux.o : common.o eflux.f90
	$(FC) $(OPT) -c eflux.f90

mhd.o : common.o mhd.f90
	$(FC) $(OPT) -c mhd.f90

urand1.o : urand1.f
	$(FC) $(OPT) -c urand1.f

clean : 
	rm *.o

cleanall : 
	rm *.o *.mod nl3d_*
