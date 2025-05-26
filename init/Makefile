################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../../..)
TESTDIR     = ../../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.launch.conf

TPROC       = 4

# required data (parameters,distributed files)
DATPARAM    = " [../../data/namelist.grads_boundary.FNL.2005053112-2015011400 namelist.grads_boundary.FNL.2005053112-2015011400] "
DATDISTS    = 


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common


all: run

