# Project: point-neuron-network-simulator

# variables for implicit rules
CXX = g++
CPPFLAGS = -std=c++11 -Wall --pedantic -Wextra -Wno-unused-parameter
CXXFLAGS = -g -O2
LDFLAGS = 
LDLIBS = -lboost_program_options
# Could faster
#   g++ -g -O2 -falign-functions=16 -falign-loops=16

# Targets and source files
BIN = bin/gen_neu
BIN_DBG = bin/gen_neu_dbg
SRCS = main.cpp math_helper.cpp legacy_lib.cpp neuron_system_utils.cpp poisson_generator.cpp
OBJS = $(SRCS:.cpp=.o)

$(BIN): $(OBJS)
	mkdir -p `dirname $(BIN)`
	$(CXX) $(LDFLAGS) -o $(BIN) $(OBJS) $(LDLIBS)

.PHONY : static-link
static-link: LDFLAGS = -static -pthread
static-link: $(BIN)

.PHONY : debug
debug: CXXFLAGS = -g
debug: $(BIN)

.PHONY : debug-O2
debug-O2: CXXFLAGS = -g -O2 -fno-omit-frame-pointer
debug-O2: $(BIN)

.PHONY : debug-log
debug-log: CPPFLAGS += -DDEBUG
debug-log: $(BIN_DBG)
$(BIN_DBG): $(BIN)
	mv $(BIN) $(BIN_DBG)

.PHONY : debug-O2-gprof
debug-O2-gprof: CXXFLAGS = -pg -g -O2
debug-O2-gprof: LDFLAGS = -pg -g -O2
debug-O2-gprof: $(BIN)

# Generate prerequisites automatically
%.d: %.cpp
	@set -e; rm -f $@; \
	  $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$

SRCDEP = $(SRCS:.cpp=.d)
include $(SRCDEP)

.PHONY : clean
clean:
	rm -f $(OBJS) $(BIN) $(SRCDEP)

