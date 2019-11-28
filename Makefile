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
MPICC = mpicc

# define compile-time flags
CFLAGS = -O3 -Wall

# define any directories containing header files
INCLUDES = -Iinc

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:
LDFLAGS =

# define any libraries to link into executable:
#   To ling libusb-1.0 :
#   LIBS = -lusb-1.0
LIBS = -lm -lopenblas
# define the source file for the library
SRC = knnring

# define the different possible executables
TYPES = sequential
MPITYPES = mpi

# define the executable file name
MAIN = main
MAINMPI = main_mpi

# define paths to .a files
LIBDIR = ./lib
# and .c files
SRCDIR = ./src

# call everytime
.PRECIOUS: %.a

all: $(addprefix $(MAIN)_, $(TYPES)) \
	$(addprefix $(MAIN)_, $(MPITYPES))

lib: $(addprefix $(LIBDIR)/, $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES)))) \
	$(addprefix $(LIBDIR)/, $(addsuffix .a, $(addprefix $(SRC)_, $(MPITYPES))))

$(MAIN)_$(TYPES): $(MAIN).c $(LIBDIR)/$(SRC)_$(TYPES).a
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

$(MAIN)_$(MPITYPES): $(MAINMPI).c $(LIBDIR)/$(SRC)_$(MPITYPES).a
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file)

$(LIBDIR)/$(SRC)_%.a: $(LIBDIR)/$(SRC)_%.o
	ar rcs $@ $<

$(LIBDIR)/$(SRC)_$(TYPES).o: $(SRCDIR)/$(SRC)_$(TYPES).c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

$(LIBDIR)/$(SRC)_$(MPITYPES).o: $(SRCDIR)/$(SRC)_$(MPITYPES).c
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	$(RM) $(LIBDIR)/*.o *~ $(addprefix $(MAIN)_, $(TYPES)) \
	$(addprefix $(MAIN)_, $(MPITYPES)) $(LIBDIR)/$(SRC)_*.a
