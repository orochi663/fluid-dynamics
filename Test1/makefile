CC=g++
CFLAGS=-c -Wl,--no-as-needed -pthread -std=c++11 `pkg-config --cflags boost` `pkg-config --cflags opencv`
LDFLAGS=-Wl,--no-as-needed -pthread `pkg-config --cflags boost` `pkg-config --libs opencv`
SOURCES=main.cpp FluidCanvas.cpp VectorField.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
