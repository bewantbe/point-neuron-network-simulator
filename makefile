# Project: vec_IFsimu

# variables for implicit rules
CXX = g++
CPPFLAGS = -std=c++11 -Wall
CXXFLAGS = -g -O2
LDFLAGS = 
LDLIBS = -lboost_program_options

# Targets and source files
BIN = bin/vec_IFsimu
SRCS = main.cpp math_helper.cpp
OBJS = $(SRCS:.cpp=.o)

$(BIN): $(OBJS)
	mkdir -p `dirname $(BIN)`
	$(CXX) $(LDFLAGS) -o $(BIN) $(OBJS) $(LDLIBS)

# Generate prerequisites automatically
%.d: %.cpp
	@set -e; rm -f $@; \
	  $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$

include $(SRCS:.cpp=.d)

.PHONY : clean
clean:
	rm $(OBJS) $(BIN)
