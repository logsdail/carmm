FC=gfortran
FCFLAGS=-c -fdefault-real-8 -O3
# LDFLAGS=-llapack
SOURCES=functions.f90 moments.f90 multipoles.f90 geometry.f90 io.f90 main.f90 
OBJECTS=$(SOURCES:.f90=.o)
EXECUTABLE=surface_dipole

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(LDFLAGS) $(OBJECTS) -o $@

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f ${OBJECTS} 

veryclean: clean
	rm -f ${EXECUTABLE}
