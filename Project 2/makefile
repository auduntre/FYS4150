CXX=g++
CXXFLAGS=-std=c++11 -larmadillo -O3
DEPS = jacobi.h
DEPSOBJ = jacobi.o
OBJ1 = buckling_beam.o
TESTS = tests-project.o
TESTMAIN = tests-main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

bb: $(DEPSOBJ) $(OBJ1)
	$(CXX) $^ $(CXXFLAGS) -o $@

tests: $(DEPSOBJ) $(TESTS) $(TESTMAIN)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	rm $(OBJ) $(DEPSOBJ) $(TESTS) tests main

cleanall:
	rm $(OBJ) $(DEPSOBJ) $(TESTS) $(TESTMAIN) tests bb
