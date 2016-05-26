CXX = g++

# CXXFLAGS = -g
CXXFLAGS = -O3

DIRSRC = ./src/
DIRFFT = $(DIRSRC)FFT/
DIRS = $(DIRSRC) $(DIRFFT)
LIBS = $(DIRSRC)source.a $(DIRFFT)fftLibrary.a
EXE = Trifid

all: $(EXE)

$(EXE): $(LIBS)
	echo Linking executable...
	$(CXX) $(CXXFLAGS) -o $(EXE) -L. $(LIBS)

$(DIRSRC)source.a: force_look
	echo Building source library...
	cd $(DIRSRC); $(MAKE) $(MFLAGS)

$(DIRFFT)fftLibrary.a: force_look
	echo Building FFT library...
	cd $(DIRFFT); $(MAKE) $(MFLAGS)

cleanObj:
	rm -f *.o *~ \#*\#
	for d in $(DIRS); do (cd $$d; $(MAKE) cleanObj); done

clean: cleanObj
	rm -f $(EXE) $(DIRSRC)*.a
	for d in $(DIRS); do (cd $$d; $(MAKE) clean); done

# Remove all images that appear to be Trifid output images.
cleanImg:
	rm -f images/*_out*

force_look:
	true
