#
# Makefile for MPI at ESIEE
#

# include local arch. particularities

LIBS    =   -lm
.c:
	mpicc -o $@ $< $(LIBS)
	mv $@ ./bin/

	
#hello: $$@.c $$@.h
#	$(CC) $(CFLAGS) -o $@ $@.c $(LDFLAGS)
#	mv $@ bin/$(PVM_ARCH)
