IDIR = ../include
CXX = g++
CXXFLAGS = -I$(IDIR) -std=c++11 -larmadillo

ODIR = obj

_DEPS = euler.h solarsystem.h celestialbody.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = euler.o solarsystem.o celestialbody.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

solarsystem: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
