# Makefile for bruits

DEPS = ./deps

lib.name = bruits
cflags += -I$(DEPS)

class.sources = ross~.c gendy~.c
gendy~.class.sources = $(DEPS)/mt19937ar/mt19937ar.c

datafiles = ross~-help.pd gendy~-help.pd

include Makefile.pdlibbuilder

.PHONY: test
test:
	$(CC) $(cflags) test_bruits.c deps/unity/unity.c -o $@
	./$@

