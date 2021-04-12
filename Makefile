CC = gcc
CPPC = g++
MPIC = mpicc
MPIPPC = mpic++

INCLUDES=-I./include -I./
CCOMPILEARG =  -Wno-unused-function -O3 -funroll-loops $(INCLUDES) -std=gnu11
CPPCOMPILEARG = -Wno-unused-function -O3 $(INCLUDES) -std=c++11

OGLLIBS = -lglfw -lGLEW -lGL -lGLU -lrt -lX11 -lXrandr -lXinerama -lXi -lXxf86vm -lXcursor
LDARGS = -lpthread -lm

SOURCES = src/NBodyHelpers.o src/NBodyInit.o src/NBodyForces.o src/NBodyKeys.o src/NBodyOctree.o src/NBodyHashedOctree.o src/NBodySerialize.o src/NBodySimulation.o src/NBodyMPISimulation.o src/NBodyMain.o
OGL_SOURCES = ogl/shader.o ogl/NBodyRenderer.o
PARALLEL_SOURCES = src/NBodyParallel.o parallel/ExecutorThreadPool.o parallel/FunctionExecutorThread.o

vpath %.cpp ./
vpath %.cpp ogl/
vpath %.o ogl/
vpath %.c ./ ./src

all : parallel 

mpi : CCOMPILEARG += -DNBODY_PARALLEL=0
mpi : CCOMPILEARG += -DNBODY_MPI=1
mpi : CPPCOMPILEARG += -DNBODY_PARALLEL=0
mpi : CPPCOMPILEARG += -DNBODY_MPI=1
mpi : CPPCOMPILEARG += -DNBODY_SIM_WITH_RENDERER=1
mpi : CC = $(MPIC)
mpi : CPPC = $(MPIPPC)
mpi : clean makeogl mpi-test.bin

mpi-noviz : CCOMPILEARG += -DNBODY_PARALLEL=0
mpi-noviz : CCOMPILEARG += -DNBODY_MPI=1
mpi-noviz : CPPCOMPILEARG += -DNBODY_PARALLEL=0
mpi-noviz : CPPCOMPILEARG += -DNBODY_MPI=1
mpi-noviz : CC = $(MPIC)
mpi-noviz : CPPC = $(MPIPPC)
mpi-noviz : clean mpi-noviz-test.bin

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

mpi-test.bin : $(SOURCES)
	$(MPIPPC) -o $@ $^ $(OGL_SOURCES) $(OGLLIBS) $(LDARGS)

mpi-noviz-test.bin : $(SOURCES)
	$(MPIPPC) -o $@ $^ $(LDARGS)

makepar : src/NBodyParallel.o 
	(cd parallel; make;)

makeogl : 
	(cd ogl; make;)

makedebugogl :
	(cd ogl; make debug;)

src/%.o: src/%.cpp include/NBodyConfig.h
	$(CPPC) -c -o $@ $(CPPCOMPILEARG) $<

src/%.o: src/%.c include/%.h
	$(CC) -c -o $@ $(CCOMPILEARG) $<


clean : 
	(cd ogl; make clean;)
	(cd parallel; make clean;)
	@rm -rf src/*.o test.bin

