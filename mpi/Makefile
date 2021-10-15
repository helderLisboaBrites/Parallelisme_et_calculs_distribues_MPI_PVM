#
# Makefile for MPI at ESIEE
#

# include local arch. particularities


.c:
	mpicc -o $@ $<
	mv $@ ./bin/

	
#hello: $$@.c $$@.h
#	$(CC) $(CFLAGS) -o $@ $@.c $(LDFLAGS)
#	mv $@ bin/$(PVM_ARCH)
