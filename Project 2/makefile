CXX=g++
CXXFLAGS=-std=c++11 -larmadillo -O3
DEPS = jacobi.h toeplitz.h
DEPSOBJ = jacobi.o toeplitz.o
OBJ1 = buckling_beam.o
OBJ2 = quantum_dots.o
TESTS = tests-project.o
TESTMAIN = tests-main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

bb: $(DEPSOBJ) $(OBJ1)
	$(CXX) $^ $(CXXFLAGS) -o $@

qd: $(DEPSOBJ) $(OBJ2)
	$(CXX) $^ $(CXXFLAGS) -o $@

tests: $(DEPSOBJ) $(TESTMAIN) $(TESTS) 
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	rm $(OBJ1) $(DEPSOBJ) $(TESTS) tests bb

cleanall:
	rm $(OBJ1) $(DEPSOBJ) $(TESTS) $(TESTMAIN) tests bb