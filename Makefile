#*
#*
#* Author: Mauricio Araya Polo
#* Date: 04/2015 - present
#*
#*

INCL=.
BIN=.

CC=gnu

#For GNU compiler
ifeq ($(CC),gnu)
	CC=gcc-9
	CFLAGS=-g -Wall -O3 -fopenmp -I$(INCL) -funroll-loops -DNOMAT -DINTEL -msse4.2
	LIBS= -lm
	#OPTS=-DGNU
endif

SRC=.

OBJS = props2D.o

%.o: $(SRC)/%.c 
	$(CC) -c $< $(CFLAGS) $(OPTS) -fpic

all: libprops

libprops: $(OBJS)
	$(CC) -o $(BIN)/$@.so $^ $(CFLAGS) $(LIBS) -shared
	#ar rcs $@.a $^
	rm -f *.o

clean:
	$(RM) *.o $(BIN)/libprops.so

cleandata:
	$(RM) *shot*.bin
