SOURCE	= mod_const.f90 mod_vars.f90 mod_fft991.f90 fft99_modified.f90 fft_d.f90 fft_i.f90 master_assim4.f90 gasdev.f90 simple_model_tl.f90 simple_model_ad.f90 simple_model.f90 hop.f90 hopt_tl.f90 hopt_ad.f90 hopL_tl.f90 hopL_ad.f90 hopt.f90
FC = gfortran
FFLAGS = -fdefault-real-8 -O3 
#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -frecursive
#FFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o main ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) *~

fort:
	 \rm -f  fort.* *~
