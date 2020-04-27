CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=cgupta_minimiser.cpp cminimiser.cpp cmorse_minimiser.cpp main.cpp Struct.cpp cXYZreadwrite.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Minimiser

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
