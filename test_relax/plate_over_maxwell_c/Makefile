# Makefile for the Plate_over_Maxwell C port
#
# all parameters are plain assignments; edit here rather than via
# environment variables

MYCC = gcc
CFLAGS = -O2 -Wall -Wextra -std=gnu99
LDFLAGS = -lm

OBJS = pom_layer.o pom_maxwell.o pom_util.o terms_2visc.o

all: pom

pom: pom_main.o $(OBJS)
	$(MYCC) -o pom pom_main.o $(OBJS) $(LDFLAGS)

pom_main.o: pom_main.c pom.h
	$(MYCC) $(CFLAGS) -c pom_main.c

pom_layer.o: pom_layer.c pom.h
	$(MYCC) $(CFLAGS) -c pom_layer.c

pom_maxwell.o: pom_maxwell.c pom.h
	$(MYCC) $(CFLAGS) -c pom_maxwell.c

pom_util.o: pom_util.c pom.h
	$(MYCC) $(CFLAGS) -c pom_util.c

# terms_2visc.c is machine generated (gen_terms.py); the file is large,
# keep optimization moderate to bound compile time
#terms_2visc.o: terms_2visc.c pom.h
#	$(MYCC) -O1 -std=gnu99 -c terms_2visc.c
terms_2visc.o: terms_2visc.c pom.h
	$(MYCC) $(CFLAGS) -c terms_2visc.c

terms_2visc.c: gen_terms.py Plate_over_Maxwell_Layer_over_Halfspace_Displacements.m
	python3 gen_terms.py

clean:
	rm -f *.o pom cmp*.txt verify_sta.txt

.PHONY: all clean
