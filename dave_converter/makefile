CC = gcc
CFLAGS = -g
ODIR = objects/$(ARCH)/
BDIR = $(HOME)/progs/bin/$(ARCH)/

SRCS = converter.c  openBC.c  openGeom.c

OBJS = $(ODIR)/converter.o $(ODIR)/openBC.o $(ODIR)/openGeom.o 

INC = converter.h

$(ODIR)/%.o:	%.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$*.o

$(BDIR)/patch2dis3d:	$(OBJS) 
	$(CC) $(CFLAGS) -o $(BDIR)/patch2dis3d $(OBJS) -lm

clean:
	/usr/bin/rm -r $(ODIR)/*.o $(ODIR)/IPA* $(ODIR)/rii*

