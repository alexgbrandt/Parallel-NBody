CC = gcc
CPPC = g++
INCLUDES=.
CCOMPILEARG = -Wall -O3 -funroll-loops -I$(INCLUDES) -std=gnu11
CPPCOMPILEARG = -Wall -O3 -I$(INCLUDES) -std=c++11

OGLLIBS = -lglfw -lGLEW -lGL -lGLU -lrt -lX11 -lXrandr -lXinerama -lXi -lXxf86vm -lXcursor
LDARGS = -lpthread -lm

SOURCES = NBodyInit.o NBodyForces.o NBodyKeys.o NBodyOctree.o NBodySimulation.o NBodyMain.o
OGL_SOURCES = ogl/shader.o ogl/NBodyRenderer.o
PARALLEL_SOURCES = NBodyParallel.o parallel/ExecutorThreadPool.o parallel/FunctionExecutorThread.o

vpath %.cpp ./
vpath %.cpp ogl/
vpath %.o ogl/


all : parallel 

serial-noviz : CCOMPILEARG += -DNBODY_PARALLEL=0
serial-noviz : CPPCOMPILEARG += -DNBODY_PARALLEL=0
serial-noviz : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=0
serial-noviz : clean noviz-serial.bin

parallel-noviz : CCOMPILEARG += -DNBODY_PARALLEL=1
parallel-noviz : CPPCOMPILEARG += -DNBODY_PARALLEL=1
parallel-noviz : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=0
parallel-noviz  : clean makepar noviz-test.bin

serial : CCOMPILEARG += -DNBODY_PARALLEL=0
serial : CPPCOMPILEARG += -DNBODY_PARALLEL=0
serial : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=1
serial : clean makeogl serial.bin

parallel : CCOMPILEARG += -DNBODY_PARALLEL=1
parallel : CPPCOMPILEARG += -DNBODY_PARALLEL=1
parallel : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=1
parallel : clean makepar makeogl test.bin

debug : CPPCOMPILEARG += -DDEBUG_OCTREE_OGL=1
debug : CCOMPILEARG += -DNBODY_PARALLEL=1
debug : CPPCOMPILEARG += -DNBODY_PARALLEL=1
debug : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=1
debug : clean makepar makedebugogl test.bin


serial.bin : $(SOURCES) 
	$(CPPC) -o test.bin $^ $(OGL_SOURCES) $(OGLLIBS) $(LDARGS)

noviz-serial.bin : $(SOURCES) 
	$(CPPC) -o test.bin $^ $(LDARGS)

noviz-test.bin: $(SOURCES)
	$(CPPC) -o test.bin $^  $(PARALLEL_SOURCES) $(LDARGS)

test.bin : $(SOURCES)
	$(CPPC) -o $@ $^ $(OGL_SOURCES) $(PARALLEL_SOURCES) $(OGLLIBS) $(LDARGS)

makepar : NBodyParallel.o 
	(cd parallel; make;)

makeogl : 
	(cd ogl; make;)

makedebugogl :
	(cd ogl; make debug;)

%.o: %.cpp NBodyConfig.h
	$(CPPC) -c $(CPPCOMPILEARG) $<

%.o: %.c %.h NBodyConfig.h
	$(CC) -c $(CCOMPILEARG) $<


clean : 
	(cd ogl; make clean;)
	(cd parallel; make clean;)
	@rm -rf *.o test.bin

