CXX=g++
CXXFLAGS=-std=c++11 -larmadillo -O3
DEPS = project.h
OBJ = project.o
TESTS = tests-project.o tests-main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

tests: $(OBJ) $(TESTS)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	rm $(OBJ) $(TESTS) tests
