CXX ?= g++
CXXFLAGS = -O3 -Wall -Wextra -lm

LDFLAGS =
INCDIRS =

DIRS = src

# Find all c and cpp files
CSRCS := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/**/*.c))
CXXSRCS := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/**/*.cpp))

# Create paths to dependency files
CDEPS := $(patsubst %.c,dep/%.c.d,$(CSRCS))
CXXEPS := $(patsubst %.cpp,dep/%.cpp.d,$(CXXSRCS))

# Create paths to object files
OBJS := \
	$(patsubst %.c,obj/%.c.o,$(CSRCS)) \
	$(patsubst %.cpp,obj/%.cpp.o,$(CXXSRCS))

all: dirs lib/libautotune.a

-include $(CDEPS)
-include $(CXXDEPS)

lib/libautotune.a: $(OBJS)
	ar rcs $@ $^

obj/%.c.o: %.c
	$(CXX) $(CXXFLAGS) -MMD -MF dep/$*.c.d -MP -c -o $@ $<

obj/%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MF dep/$*.cpp.d -MP -c -o $@ $<

clean:
	-rm -r dep/* obj/* lib/*

dirs:
	mkdir -p $(foreach DIR,$(DIRS),lib \
		$(patsubst %,dep/%,$(wildcard $(DIR)/**)) \
		$(patsubst %,obj/%,$(wildcard $(DIR)/**)))


.PHONY: dirs clean
