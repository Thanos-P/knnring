# ####################################################################
#
#			   C/C++ Makefile
#
#	Modified by: Thanos Paraskevas
#							 athanasps <athanasps@ece.auth.gr>
# Original:
# Author: Dimitris Floros <fcdimitr@auth.gr>
#
# Adapted from
#  http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#
# ####################################################################
#
# 'make'        build executable file 'main'
# 'make lib'	build the libraries .a
# 'make clean'  removes all .o and executable files
#

# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7

# define compile-time flags
CFLAGS = -O3 -Wall

# define any directories containing header files
INCLUDES =

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:
LDFLAGS =

# define any libraries to link into executable:
#   To ling libusb-1.0 :
#   LIBS = -lusb-1.0
LIBS = -lm

# define the source file for the library
SRC = kNN

# define the different possible executables
TYPES = sequential

# define the executable file  name
MAIN = main

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above
#

# call everytime
.PRECIOUS: %.a

all: $(addprefix $(MAIN)_, $(TYPES))

lib: $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES)))

$(MAIN)_%: $(MAIN).c $(SRC)_%.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file)

.o.a:
	ar rcs $@ $<

# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	$(RM) *.o *~ $(addprefix $(MAIN)_, $(TYPES)) *.a
