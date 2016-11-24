CC        = gcc
CFLAGS   += -std=c99
CFLAGS   += -Wall -O3
CFLAGS   += -DT_MIN=5
CFLAGS   += -DCOLMAJOR=1
CFLAGS   += -DUPPER=1
CFLAGS   += -DUNIT=0
CFLAGS   += -DDGEMV_DOTF_FUSE=5
CFLAGS   += -DDGEMV_AXPYF_FUSE=4

Arch := intel64
IntelDir := /opt/intel/compilers_and_libraries/linux
IncludeDir := $(IntelDir)/mkl/include
LibDirs := $(IntelDir)/lib/$(Arch) $(IntelDir)/mkl/lib/$(Arch)
LinkerFlag := -Wl,
Rpath := $(patsubst %,$(LinkerFlag)-rpath $(LinkerFlag)%,$(LibDirs))
Lpath := $(patsubst %,-L%,$(LibDirs))

LDFLAGS += $(Lpath) $(Rpath)
CPPFLAGS += -I$(IncludeDir) -DMKL_ILP64
LDLIBS += -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lm -lpthread

TARGETS = gemv

all: $(TARGETS)

clean:
	$(RM) $(TARGETS)

plot:
	gnuplot gemv_time.gnuplot
	gnuplot gemv_MFLOPS.gnuplot

bench:
	./gemv > gemv.data
