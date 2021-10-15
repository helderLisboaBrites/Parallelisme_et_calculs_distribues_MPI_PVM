#
# minimum Makefile for PVM at ESIEE
#

# include local arch. particularities

include $(PVM_ROOT)/conf/$(PVM_ARCH).def

CFLAGS	=	-I$(PVM_ROOT)/include $(ARCHCFLAGS)

LIBS	=	-lpvm3 -lgpvm3 $(ARCHLIB) -lm
LFLAGS	=	-L$(PVM_ROOT)/lib/$(PVM_ARCH)
LDFLAGS	=	$(LFLAGS) $(LIBS)

# dependances go here

.c:
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
	mv $@ $(PVM_ROOT)/bin/$(PVM_ARCH)

#hello: $$@.c $$@.h
#	$(CC) $(CFLAGS) -o $@ $@.c $(LDFLAGS)
#	mv $@ bin/$(PVM_ARCH)
