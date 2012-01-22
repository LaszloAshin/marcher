CPPFLAGS := -march=native -pipe -fomit-frame-pointer -Wall -pedantic -pthread
CXXFLAGS := -Weffc++

marcher: CPPFLAGS := $(CPPFLAGS) -O3

shader: LDLIBS := -lSDL -lGL
shader: CPPFLAGS := $(CPPFLAGS) -Os -s

BINS := marcher shader

all: $(BINS)

clean:
	rm -f $(BINS)
