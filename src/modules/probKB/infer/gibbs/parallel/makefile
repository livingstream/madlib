OBJS = gibbs_scheduler.o gibbs_gist.o gibbs_cgla.o test_gibbs_gist.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

test : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o test -lpthread 

gibbs_scheduler.o : gibbs_scheduler.h gibbs_scheduler.cpp
	$(CC) $(CFLAGS) gibbs_scheduler.cpp

gibbs_gist.o : gibbs_gist.h gibbs_gist.cpp
	$(CC) $(CFLAGS) gibbs_gist.cpp

gibbs_cgla.o : gibbs_cgla.h gibbs_cgla.cpp
	$(CC) $(CFLAGS) gibbs_cgla.cpp

test_gibbs_gist.o : test_gibbs_gist.cpp
	$(CC) $(CFLAGS) test_gibbs_gist.cpp
clean: 
	rm *.o test
