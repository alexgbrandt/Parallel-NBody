CC = gcc
CPPC = g++
INCLUDES=.
CCOMPILEARG = -Wall -g -I$(INCLUDES) -std=gnu11
CPPCOMPILEARG = -Wall -g -I$(INCLUDES) -std=c++11
LDARGS = -lglfw -lGLEW -lGL -l GLU -lrt -lX11 -lXrandr -lXinerama -lXi -lXxf86vm -lXcursor -lm

SOURCES = shader.o NBodyRenderer.o

all : $(SOURCES)

debug : CPPCOMPILEARG += -DDEBUG_OCTREE_OGL=1
debug : all

%.o: %.cpp %.hpp
	$(CPPC) -c $(CPPCOMPILEARG) $<

%.o: %.c %.h
	$(CC) -c $(CCOMPILEARG) $<

clean :
	@rm -rf *.o ogl-test.bin


