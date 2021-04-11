CC     = g++
CFLAGS = -g -w -Iinc/
LIBS   = #-lSDL2 -lSDL2_image -lSDL2_ttf

SRCS   = $(filter-out src/image_drawing.cpp, $(wildcard src/*.cpp))
OBJS   = $(patsubst src/%.cpp,bin/%.o,$(SRCS))
DEPS   = $(OBJS:.o:=.d)
DIRS   = src inc bin
EXE    = a.out

all: $(DIRS) $(EXE)

$(DIRS):
	mkdir -p $@

$(EXE): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

bin/%.o : src/%.cpp
	$(CC) -o $@ $(CFLAGS) -c $<

bin/%.o : src/%.cpp inc/%.h
	$(CC) -o $@ $(CFLAGS) -c $<

run : all
	./$(EXE)

clean:
	rm -rf bin *~ $(EXE)
