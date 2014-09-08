CC=g++
CFLAGS=-c -Wl,--no-as-needed -pthread -std=c++11 `pkg-config --cflags boost` `pkg-config --cflags opencv`
LDFLAGS=-Wl,--no-as-needed -pthread `pkg-config --cflags boost` `pkg-config --libs opencv`
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) -o $@ $^

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)
