all: release

release: Atomaker
debug: Atomaker

%.d: %.cpp makefile
	$(CXX) $(CXXFLAGS) -MM -o $@ $<

CXXFLAGS += -std=c++1y
release: CXXFLAGS += -O3
debug: CXXFLAGS += -g3

-include main.d
-include main_d.d


OBJECTS = main.o
Atomaker: $(OBJECTS) $(LDLIBS)
	$(CXX) $(LDFLAGS) -o $@ $<
clean:
	$(RM) Atomaker $(OBJECTS)
.PHONY: clean all release debug
.SECONDARY: $(OBJECTS:%.o=%.d %.o)
