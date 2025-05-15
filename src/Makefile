SOURCE	= mod_const.f90 mod_vars.f90 mod_fft991.f90 fft99_modified.f90 fft_d.f90 fft_i.f90 master_assim.f90 gasdev.f90 burgers2_tl.f90 burgers2_ad.f90 burgers2.f90 hopt.f90 hopt_tl.f90 hopt_ad.f90 chavarin.f90 chavarin_ad.f90 cost_function.f90 init.f90 create_obs.f90 congrad.f90 correl.f90
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
