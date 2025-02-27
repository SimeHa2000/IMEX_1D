# Compiler
CXX = clang++

# Flags
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra

# Source and Header files
SRCS = main.cpp parms.cpp IMEX.cpp EoS.cpp
HEADERS = parms.h IMEX.h EoS.h
OBJS = main.o parms.o IMEX.o EoS.o

# Target: Compile and link everything
IMEX_1d: $(OBJS)
	$(CXX) $(CXXFLAGS) -o IMEX_1d $(OBJS)

# Rule for compiling object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS)
