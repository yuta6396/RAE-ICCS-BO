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
PPCONF      = pp.d01.conf

TPROC       = 4

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common


all: run

