# Basic makefile for tempo2 plugins to be compiled outside of tempo2

# Get system arch automatically
ARCH = $(shell uname | tr A-Z a-z)
# Or just hard-code here if you're having trouble, eg:
#ARCH = linux

# This is the path to your copy of the tempo2 source.
# Note, in general different from your $TEMPO2 directory.
TEMPO2_SRC = /data1/demorest/src/tempo2

# Include dirs, etc
CXXFLAGS = -O3 -I $(TEMPO2)/include -I $(TEMPO2_SRC)

# Extra libraries
LIBS = -L $(TEMPO2)/lib -ltempo2 -lcfitsio -lm 

# Base name of all plugins to be compiled
PLUGINS = mcmc grid

PLUGIN_SRCS = $(addsuffix _plug.C, $(PLUGINS))
PLUGIN_LIBS = $(addsuffix _$(ARCH)_plug.t2, $(PLUGINS))

all: $(PLUGIN_LIBS)

install: $(PLUGIN_LIBS)
	cp -p $(PLUGIN_LIBS) $(TEMPO2)/plugins

clean:
	rm -f $(PLUGIN_LIBS)

%_$(ARCH)_plug.t2: %_plug.C
	$(CXX) $(CXXFLAGS) -fPIC -shared -o $@ $< $(LIBS)

