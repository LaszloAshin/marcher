CPPFLAGS := -march=native -O3 -pipe -fomit-frame-pointer -Wall -pedantic -Weffc++ -pthread

BINS := marcher

all: $(BINS)

clean:
	rm -f $(BINS)
