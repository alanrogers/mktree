CFLAGS = -Wall -g
#CFLAGS = -ggdb -Wall -Wshadow -Wstrict-prototypes
CC = gcc
.c.o:
	$(CC) $(CFLAGS) -c $<
TARGETS = mktree mkini.py
BINDIR = $(HOME)/bin

MKTREE = mktree.o textree.o iscoales.o unirand.o memstack.o bye.o mutation.o \
      getcic.o alloc2d.o
mktree : $(MKTREE) Makefile
	$(CC) $(CFLAGS) -o $@ $(MKTREE) -lm

BYHAND = byhand.o textree.o iscoales.o unirand.o memstack.o bye.o \
  mutation.o getcic.o alloc2d.o
byhand : $(BYHAND) Makefile
	$(CC) $(CFLAGS) -o $@ $(BYHAND) -lm

clean :
	rm -f *.o $(TARGET) depend *.aux *.dvi *.log *~

mktree.o : mktree.c iscoales.h Makefile memstack.h unirand.h mytypes.h


byhand.o : byhand.c iscoales.h Makefile memstack.h mytypes.h

memstack.o : memstack.c memstack.h Makefile

unirand.o : unirand.c unirand.h mytypes.h Makefile

install : $(TARGETS)
	cp $(TARGETS) $(BINDIR)
	touch install

depend :
	echo '#Automatically generated dependency info' > depend
	$(CC) -M *.c >> depend

include depend
