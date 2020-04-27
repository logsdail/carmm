CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp fcc_neighbours.cpp calculations.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=polynomial

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
