# Makefile for bruits


lib.name = bruits
cflags += -I./deps

class.sources = ross~.c gendy~.c
gendy~.class.sources = deps/mt19937ar/mt19937ar.c

datafiles = ross~-help.pd gendy~-help.pd

include Makefile.pdlibbuilder
