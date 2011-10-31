CPPFLAGS := -march=native -O3 -pipe -fomit-frame-pointer

BINS := marcher

all: $(BINS)

clean:
	rm -f $(BINS)
