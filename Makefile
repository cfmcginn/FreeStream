CXX=g++
#O3 for max optimization (go to 0 for debug)
#Turn off -Werror
CXXFLAGS=-Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS+=-Wno-error=misleading-indentation
endif

INCLUDE=-I$(PWD)

SPHENIXPATH=/opt/sphenix/core/gsl/include

ifneq ("$(wildcard $(SPHENIXPATH))","")
  $(echo ADDING $(SPHENIXPATH))
  INCLUDE+=-I/opt/sphenix/core/gsl/include -L/opt/sphenix/core/gsl/lib/
endif

FDASH=-fpermissive -fopenmp
LDASH=-lgsl -lgslcblas

ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_LIB=mkdir -p $(PWD)/lib
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir
MKDIR_LOG=mkdir -p $(PWD)/logdir

all: mkdirBin mkdirLib mkdirPdf mkdirlogdir mkdirOutput lib/paramreader.o lib/linear_int.o bin/initE.exe bin/FS.exe bin/InitED.exe bin/initPointSource.exe bin/processED.exe bin/freestream.exe bin/initIPGlasma.exe bin/compRootToDat.exe bin/compareIPFS.exe bin/compDat.exe bin/checkInitFiles.exe bin/rescaleInitED.exe bin/checkED.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

mkdirlogdir:
	$(MKDIR_LOG)

lib/paramreader.o: src/paramreader.cpp
	$(CXX) $(CXXFLAGS) -c src/paramreader.cpp $(INCLUDE) -o lib/paramreader.o

lib/linear_int.o: src/linear_int.cpp
	$(CXX) $(CXXFLAGS) -c src/linear_int.cpp $(INCLUDE) -o lib/linear_int.o

bin/initE.exe: src/initE.cpp
	$(CXX) $(CXXFLAGS) $(FDASH) src/initE.cpp $(LDASH) lib/paramreader.o $(INCLUDE) -o bin/initE.exe

bin/FS.exe: src/FS.cpp
	$(CXX) $(CXXFLAGS) $(FDASH) src/FS.cpp $(LDASH) lib/paramreader.o lib/linear_int.o $(INCLUDE) -o bin/FS.exe

bin/InitED.exe: src/InitED.C
	$(CXX) $(CXXFLAGS) src/InitED.C $(INCLUDE) $(ROOT) -o bin/InitED.exe

bin/initPointSource.exe: src/initPointSource.C
	$(CXX) $(CXXFLAGS) src/initPointSource.C $(INCLUDE) $(ROOT) -o bin/initPointSource.exe

bin/initIPGlasma.exe: src/initIPGlasma.C
	$(CXX) $(CXXFLAGS) src/initIPGlasma.C $(INCLUDE) $(ROOT) -o bin/initIPGlasma.exe

bin/processED.exe: src/processED.C
	$(CXX) $(CXXFLAGS) src/processED.C $(INCLUDE) $(ROOT) -o bin/processED.exe

bin/freestream.exe: src/freestream.C
	$(CXX) $(CXXFLAGS) src/freestream.C $(INCLUDE) $(ROOT) -o bin/freestream.exe

bin/compRootToDat.exe: src/compRootToDat.C
	$(CXX) $(CXXFLAGS) src/compRootToDat.C $(INCLUDE) $(ROOT) -o bin/compRootToDat.exe

bin/compareIPFS.exe: src/compareIPFS.C
	$(CXX) $(CXXFLAGS) src/compareIPFS.C $(INCLUDE) $(ROOT) -o bin/compareIPFS.exe

bin/compDat.exe: src/compDat.C
	$(CXX) $(CXXFLAGS) src/compDat.C $(INCLUDE) $(ROOT) -o bin/compDat.exe

bin/checkInitFiles.exe: src/checkInitFiles.C
	$(CXX) $(CXXFLAGS) src/checkInitFiles.C $(INCLUDE) $(ROOT) -o bin/checkInitFiles.exe

bin/rescaleInitED.exe: src/rescaleInitED.C
	$(CXX) $(CXXFLAGS) src/rescaleInitED.C $(INCLUDE) $(ROOT) -o bin/rescaleInitED.exe

bin/checkED.exe: src/checkED.C
	$(CXX) $(CXXFLAGS) src/checkED.C $(INCLUDE) -o bin/checkED.exe

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f lib/*.o
	rm -rf lib
	rm -f configs/*~
	rm -f configs/#*#
	rm -f data/*~
	rm -f data/#*#
	rm -f include/*~
	rm -f include/#*#
	rm -f input/*~
	rm -f input/#*#
	rm -f src/*~
	rm -f src/#*#
