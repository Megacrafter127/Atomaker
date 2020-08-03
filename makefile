all: Atomaker

%.d: %.cpp makefile
	$(CXX) $(CXXFLAGS) -MM -o $@ $<

-include main.d

CXXFLAGS += -O3 -std=c++1y

OBJECTS = main.o

Atomaker: $(OBJECTS) $(LDLIBS)
	$(CXX) $(LDFLAGS) -o $@ $<
clean:
	$(RM) $(OBJECTS)
.PHONY: clean all
.SECONDARY: $(OBJECTS:%.o=%.d %.o)
