CC = gcc
CPPC = g++
INCLUDES=./include
CCOMPILEARG = -Wall -O3 -funroll-loops -I$(INCLUDES) -std=gnu11
CPPCOMPILEARG = -Wall -O3 -funroll-loops -I$(INCLUDES) -std=c++11

SOURCES = ExecutorThreadPool.o FunctionExecutorThread.o

vpath %.cpp ./src 

all : $(SOURCES)

%.o: %.cpp 
	$(CPPC) -c $(CPPCOMPILEARG) $<

%.o: %.c %.h
	$(CC) -c $(CCOMPILEARG) $<

clean : 
	@rm -rf *.o 


