CXX = g++
CXXFLAGS = -O3
OUTPUT = fftLibrary.a
INCLUDE = -I . -I ..

OBJS = dxpose.o fft2d.o fftext.o fftlib.o matlib.o

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
