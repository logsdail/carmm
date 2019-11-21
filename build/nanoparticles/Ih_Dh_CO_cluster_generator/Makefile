CC=g++
CFLAGS=-c -Wall -pedantic
LDFLAGS=
SOURCES=COct.cpp CCub.cpp CDec.cpp CGeometry.cpp CIco.cpp CPentBi.cpp Struct.cpp cXYZreadwrite.cpp geometrymaker.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=structure_gen

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
