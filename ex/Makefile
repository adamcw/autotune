# Blossom V
BLOSSOMV = ../blossomv
BDIRS := $(BLOSSOMV) $(BLOSSOMV)/MinCost $(BLOSSOMV)/GEOM
BSOURCES := $(foreach dir, $(BDIRS), $(wildcard $(dir)/*.cpp))
# Filter out example.cpp as a source
BSOURCES := $(filter-out $(BLOSSOMV)/example.cpp, $(BSOURCES))
BOBJS := $(patsubst %.cpp, %.o, $(BSOURCES))

# AutoTune
LIBAUTOTUNE = ../libautotune

LDFLAGS = -L $(LIBAUTOTUNE)/lib -l autotune
INCDIRS = -I $(LIBAUTOTUNE)/src

# Add librt if the target platform is Linux
ifeq ($(shell uname -s),Linux)
    LDFLAGS += -lrt
endif

CXX ?= g++

ifeq ($(ATDEBUG),1)
CXXFLAGS = -g3 -Wall -Wextra -lm $(INCDIRS)
else
CXXFLAGS = -O3 -Wall -Wextra -lm $(INCDIRS)
endif

DIRS = src

# Find all c and cpp files
CSRCS := $(foreach DIR, $(DIRS), $(wildcard $(DIR)/*.c))
CXXSRCS := $(foreach DIR, $(DIRS), $(wildcard $(DIR)/*.cpp))

# Create paths to dependency files
CDEPS := $(patsubst %.c, dep/%.c.d, $(CSRCS))
CXXDEPS := $(patsubst %.cpp, dep/%.cpp.d, $(CXXSRCS))

# Create paths to object files
OBJS := \
	$(patsubst %.c, obj/%.c.o, $(CSRCS)) \
	$(patsubst %.cpp, obj/%.cpp.o, $(CXXSRCS))

# Binaries 
BINS := \
	$(foreach DIR, $(DIRS), $(patsubst $(DIR)/%.c, bin/%, $(CSRCS))) \
	$(foreach DIR, $(DIRS), $(patsubst $(DIR)/%.cpp, bin/%, $(CXXSRCS)))

# Build all files in the /src/ directory
all: dirs $(BINS)

-include $(CDEPS)
-include $(CXXDEPS)

# Build a specific binary 
bin/%: obj/src/%.c.o $(LIBAUTOTUNE)/lib/libautotune.a 
	$(CXX) -o $@ $^ $(BOBJS) $(LDFLAGS)

# Make a c object file
obj/%.c.o: %.c
	$(CXX) $(CXXFLAGS) -MMD -MF dep/$*.c.d -MP -c -o $@ $<

# Make a cpp object file
obj/%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MF dep/$*.cpp.d -MP -c -o $@ $< 

# Remove all compiled content
clean:
	-rm -r bin/* dep/* obj/*; 

# Build all required folders
dirs:
	mkdir -p $(foreach DIR,$(DIRS),bin dep/$(DIR) obj/$(DIR))

.PHONY: dirs clean 

