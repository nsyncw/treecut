CCOMP = gcc
CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm
#CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm

all: hi_treem 
hi_treem: hi_treem.c parser_treem.c timer.c
	$(CCOMP) $(CFLAGS) -o hi_treem hi_treem.c 
clean: 
	rm -f hi_treem *~
