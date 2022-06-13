# SpM Makefile
include ./make.inc

.SUFFIXES : .mod .o .f90

%.o : %.f90
	$(FC) $(FCFLAGS) -c $< $(LIB) 

			
SRC =	lib_math.f90 	\
			spm_util.f90	\
			spm_svd.f90		\
			spm_admm.f90	\
			spm_main.f90

OBJ = $(SRC:.f90=.o)


default : xspm.x

xspm.x : $(OBJ)
	$(FC) $(FCFLAGS) $(OBJ) -o $@ $(LIB)

clean :
		rm -f *.o *.mod *~

veryclean :
		rm -f *.o *.mod *~ *.x 
