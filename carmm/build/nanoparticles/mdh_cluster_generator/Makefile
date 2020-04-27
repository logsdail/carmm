FC=gfortran
FCFLAGS=-c -fdefault-real-8 -O3
# LDFLAGS=-llapack
SOURCES=marks3.for 
OBJECTS=$(SOURCES:.f90=.o)
EXECUTABLE=mdh_generator

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
