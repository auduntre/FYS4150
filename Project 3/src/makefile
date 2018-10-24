IDIR = ../include
CXX = g++
CXXFLAGS = -I$(IDIR) -std=c++11 -larmadillo -O3 #-ffast-math #-fopt-info-all

ODIR = obj

_DEPS = euler.h verlet.h solarsystem.h celestialbody.h betasolarsystem.h \
		gensolarsystem.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = euler.o verlet.o solarsystem.o celestialbody.o betasolarsystem.o \
	   gensolarsystem.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_TBOBJ = twobody.o
TBOBJ = $(patsubst %,$(ODIR)/%,$(_TBOBJ))

_TBMOBJ = twobodymercury.o
TBMOBJ = $(patsubst %,$(ODIR)/%,$(_TBMOBJ))

_TBBOBJ = twobodybeta.o
TBBOBJ = $(patsubst %,$(ODIR)/%,$(_TBBOBJ))

_TME = timemeasure.o
TME = $(patsubst %,$(ODIR)/%,$(_TME))



$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


timemeasure: $(OBJ) $(TME)
	$(CXX) -o $@ $^ $(CXXFLAGS)


twobody: $(OBJ) $(TBOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)


twobodymercury: $(OBJ) $(TBMOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)


twobodybeta: $(OBJ) $(TBBOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~
