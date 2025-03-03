# Compiler
CXX = clang++

#Libraries
Eigen_DIR = /Users/simeon/Desktop/physics_phd/code/bin/eigen

# Flags
CXXFLAGS = -std=c++14 -O3 -Wall -Wextra -I$(Eigen_DIR)

# Source and Header files
SRCS = main.cpp parms.cpp IMEX.cpp EoS.cpp BoundaryConditions.cpp
HEADERS = parms.h IMEX.h EoS.h BoundaryConditions.h
OBJS = main.o parms.o IMEX.o EoS.o BoundaryConditions.o

# Target: Compile and link everything
IMEX_1d: $(OBJS)
	$(CXX) $(CXXFLAGS) -o IMEX_1d $(OBJS)

# Rule for compiling object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS)