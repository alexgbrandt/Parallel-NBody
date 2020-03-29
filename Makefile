CC = gcc
CPPC = g++
INCLUDES=.
CCOMPILEARG = -Wall -g -funroll-loops -I$(INCLUDES) -std=gnu11
CPPCOMPILEARG = -Wall -g -I$(INCLUDES) -std=c++11
LDARGS = -lglfw -lGLEW -lGL -l GLU -lrt -lX11 -lXrandr -lXinerama -lXi -lXxf86vm -lXcursor -lpthread -lm

SOURCES = NBodyInit.o NBodyForces.o NBodyKeys.o NBodyOctree.o NBodyParallel.o NBodySimulation.o NBodyMain.o
OGL_SOURCES = ogl/shader.o ogl/NBodyRenderer.o
PARALLEL_SOURCES = parallel/ExecutorThreadPool.o parallel/FunctionExecutorThread.o

vpath %.cpp ./ 
vpath %.cpp ogl/
vpath %.o ogl/


all : makepar makeogl test.bin

debug : CPPCOMPILEARG += -DDEBUG_OCTREE_OGL=1
debug : makedebugogl test.bin

makepar : 
	(cd parallel; make;)

makeogl : 
	(cd ogl; make;)

makedebugogl :
	(cd ogl; make debug;)

%.o: %.cpp NBodyConfig.h
	$(CPPC) -c $(CPPCOMPILEARG) $<

%.o: %.c %.h NBodyConfig.h
	$(CC) -c $(CCOMPILEARG) $<

test.bin : $(SOURCES) 
	$(CPPC) -o $@ $^ $(OGL_SOURCES) $(PARALLEL_SOURCES) $(LDARGS)

clean : 
	(cd ogl; make clean;)
	(cd parallel; make clean;)
	@rm -rf *.o test.bin


