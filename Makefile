# MAKEFILE FOR simple C++ programming

CFLAGS = -O2 -g -pedantic -Wall -std=c++11
INCLUDE = -I./src -I./src/Departure -I./include/
LIBS = -L./lib/
CXX = g++
SRCDIR = src
HDRDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp
HDREXT = hpp
OBJEXT = o

SUBDIRS = Departure

$(shell mkdir -p $(OBJDIR))
$(shell mkdir -p $(addprefix $(OBJDIR)/, $(SUBDIRS)))

FILES = main \
Departure/DepartureCoeffs
				
SRCS = $(FILES:=.$(SRCEXT))
HDRS = $(FILES:=.$(HDREXT))
OBJS = $(FILES:=.$(OBJEXT))
#HDRS = $(SRCS:.$(SRCEXT)=.$(HDREXT))
#OBJS = $(SRCS:.$(SRCEXT)=.$(OBJEXT))
FULLPATHOBJ = $(addprefix $(OBJDIR)/, $(OBJS))
FULLPATHSRC = $(addprefix $(SRCDIR)/, $(SRCS))
FULLPATHHDR = $(addprefix $(HDRDIR)/, $(HDRS))

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(HDRDIR)/%.$(HDREXT)
	mkdir -p $(dir $@)
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)
	
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	mkdir -p $(dir $@)
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)

rrldepart : $(FULLPATHOBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBS)

.PHONY: clean
clean:
	rm $(FULLPATHOBJ) rrldepart
