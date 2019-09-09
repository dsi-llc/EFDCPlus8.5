.SUFFIXES:

.SUFFIXES: .f90 .F90 .for .FOR .o


INCLUD   =  . 
CDFDIR   =  .
CDFLIBS  =  

### LINUX Intel Compiler call ############################################

FC       = /opt/intel/bin/ifort 

FFLAGS   = -O3
#FFLAGS  = -O
#FFLAGS  = -g  
#FFLAGS  = -g -C
#FFLAGS  = -fast

########################################################################

PROGRAM = getefdc
PROGEXT = .x

PROJHOME= ..

DATADIR = $(PROJHOME)/SHED


# PROGRAM COMPONENTS  ##################################################

OBJS    =  infomod.o  efdcpromod.o tecmod.o geteeoutmod.o       \
           xyijconv.o  gethfreqout.o 

MODS    =  globalvars.o 

MAIN    =  getefdc.o

EXEPROG = $(PROGRAM)$(PROGEXT)

###### FOR LINUX #######################################################	

$(EXEPROG) : $(OBJS) $(MODS) $(MAIN) 
  
	$(FC) $(FFLAGS) $(CDFDIR) $(CDFLIBS) $(INCLUD) -o  $@  $(OBJS) $(MODS) $(MAIN) 

# Recompile only source file of change #################################

.f90.o: 
	$(FC) $(FFLAGS)  $(INCLUD) -c $<

.F90.o: 
	$(FC) $(FFLAGS)  $(INCLUD) -c $<

.for.o: 
	$(FC) $(FFLAGS)  $(INCLUD) -c $<

.FOR.o: 
	$(FC) $(FFLAGS)  $(INCLUD) -c $<

### Special dependencies ###############################################

$(OBJS)              : $(MODS) 
$(MAIN)              : $(OBJS)
 
clean:
	rm -f *.o *.mod *.x
