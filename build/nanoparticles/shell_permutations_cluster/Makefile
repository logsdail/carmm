CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp cPermute.cpp Struct.cpp cXYZreadwrite.cpp cgupta_minimiser.cpp cminimiser.cpp cmorse_minimiser.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Permute

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
