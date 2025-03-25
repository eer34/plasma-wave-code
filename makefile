########################################################################
#     Makefile
########################################################################
#module add Intel_compiler/19.0.5
#module add MPI/Intel/IMPI
#module add MKL

FC          = mpiifort
MAIN        = IWNS
# FCC         = -ffree-line-length-none
FCC         =
INCLUDE     = 
LIB         =
# INCLUDE     = -I/public/software/compiler/intel-compiler/2021.3.0/mkl/include
# LIB         =-L/public/software/compiler/intel-compiler/2021.3.0/mkl/lib/intel64
FCCFLAG     = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
GRAPH_PATH  = ./eps/
OBJECT_PATH = ./obj/
TEXT_PATH   = ./txt/

OBJECTS = errfun.o m_spline_interp.o quadpack_generic.o mathfunction.o tool.o   zpl5.o IWNS.o

# Makefile
${MAIN} : ${OBJECTS}
	  ${FC} ${FCC} -o ${MAIN} ${OBJECTS} ${INCLUDE}${LIB} ${FCCFLAG}
errfun.o :errfun.f90
	   ${FC} ${FCC} -c errfun.f90
m_spline_interp.o   : m_spline_interp.f90
	   ${FC} ${FCC} -c m_spline_interp.f90
quadpack_generic.o : quadpack_generic.F90
	   ${FC} ${FCC} -c quadpack_generic.F90
mathfunction.o: errfun.f90 quadpack_generic.F90 mathfunction.f90
	   ${FC} ${FCC} -c mathfunction.f90 -mkl
tool.o   :errfun.f90 mathfunction.f90 tool.f90
	   ${FC} ${FCC} -c tool.f90 -mkl
zpl5.o   : m_spline_interp.f90 quadpack_generic.F90 zpl5.f90
	   ${FC} ${FCC} -c zpl5.f90 -mkl
IWNS.o   :errfun.f90 mathfunction.f90 tool.f90 m_spline_interp.f90 quadpack_generic.F90 zpl5.f90 IWNS.f90
	   ${FC} ${FCC} -c IWNS.f90 -mkl
move :
	-mv -f $(MAIN) *.o *.mod ${OBJECT_PATH}
	-mv -f *.txt ${TEXT_PATH}
	-mv -f *.eps ${GRAPH_PATH}

clean :
	-rm -f ${OBJECT_PATH}*
	-rm -f ${TEXT_PATH}*
	-rm -f *.exe
	-rm -f *.err
	-rm -f *.out
	-rm -f *.o
	-rm -f *.mod
clear:
	-rm -f *.txt