CXX = g++

# CXXFLAGS = -g
CXXFLAGS = -O3

OUTPUT = source.a
INCLUDE = -I . -I ./FFT

OBJS = ImageInterface.o ImageReal.o Image3ChnlReal.o Random.o Bayer.o DebayerProcs.o ImageProcs.o FourierTransform.o SteerablePyramid.o ImagePairStats.o main.o

all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	ar rv $(OUTPUT) $(OBJS)
	ranlib $(OUTPUT)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

cleanObj:
	rm -f *.o *~ \#*\#

clean: cleanObj
	rm -f *.a

force_look:
	true
