appname := app

CXX := clang++
CXXFLAGS := -std=c++11 -O3

PARALLEL_RUN=${PARALLEL}

ifeq ($(PARALLEL_RUN),)
# takes values 0 or 1
PARALLEL_RUN="1"
endif

MYFLAGS = -DPARALLEL_RUN=$(PARALLEL_RUN)
CXXFLAGS += $(MYFLAGS)

ifeq ($(PARALLEL_RUN),"1")
# Link in the TBB libraries
LDLIBS := -ltbb
else
LDLIBS :=
endif

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend