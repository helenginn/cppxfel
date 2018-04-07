BOOST_HEADER_PATH=/well/strubi/ginn/cctbx2/modules/boost/
CCP4_HEADER_PATH=/apps/strubi/ccp4/ccp4-6.5/include/
HDF5_HEADER_PATH=/well/strubi/ginn/cctbx2/base/include/
PNG_HEADER_PATH=/well/strubi/ginn/cctbx2/base/include/

BOOST_LIBRARY_PATH=/well/strubi/ginn/cctbx2/build_tmp/lib/
CCP4_LIBRARY_PATH=/apps/strubi/ccp4/ccp4-6.5/lib/
HDF5_LIBRARY_PATH=/well/strubi/ginn/cctbx2/base/lib/
PNG_LIBRARY_PATH=/well/strubi/ginn/cctbx2/base/lib/

CPPXFEL_HEADER_PATHS=-I $(PNG_HEADER_PATH) -I $(BOOST_HEADER_PATH) -I $(CCP4_HEADER_PATH) -I $(HDF5_HEADER_PATH)

CPPXFEL_LIBRARY_PATHS=-L $(PNG_HEADER_PATH) -L $(BOOST_LIBRARY_PATH) -L $(CCP4_LIBRARY_PATH) -L $(HDF5_LIBRARY_PATH)

FLAGS = -O3 -std=c++0x

-rpath=/well/strubi/ginn/cctbx2/base/lib/

#-Wl,-rpath=/well/strubi/ginn/cctbx2/base/lib/

LIBFLAGS = -lm -lpng12 -l:libhdf5.so.10 -lhdf5_hl -lccp4c -lboost_system-mt -lboost_thread-mt

#LDLIBS=-Wl,-rpath,/well/strubi/ginn/cctbx2/base/lib/ -lpng
#LDFLAGS = -static

BEFORE = $(CPPXFEL_HEADER_PATHS) $(FLAGS)
AFTER = $(CPPXFEL_LIBRARY_PATHS) $(LIBFLAGS)

all:
	@echo cppxfel header paths = $(CPPXFEL_HEADER_PATHS)
	@echo cppxfel library paths = $(CPPXFEL_LIBRARY_PATHS)
	@echo flags = $(FLAGS)
	@echo library flags = $(LIBFLAGS)

	g++ $(BEFORE) -c AmbiguityBreaker.cpp
	g++ $(BEFORE) -c Beam.cpp
	g++ $(BEFORE) -c CSV.cpp
	g++ $(BEFORE) -c Detector.cpp
	g++ $(BEFORE) -c FileParser.cpp
	g++ $(BEFORE) -c FileReader.cpp
	g++ $(BEFORE) -c FreeLattice.cpp
	g++ $(BEFORE) -c GaussianBeam.cpp
	g++ $(BEFORE) -c GeometryParser.cpp
	g++ $(BEFORE) -c GeometryRefiner.cpp
	g++ $(BEFORE) -c GraphDrawer.cpp
	g++ $(BEFORE) -c Hdf5Crystal.cpp
	g++ $(BEFORE) -c Hdf5Image.cpp
	g++ $(BEFORE) -c Hdf5Manager.cpp
	g++ $(BEFORE) -c Hdf5ManagerCheetah.cpp
	g++ $(BEFORE) -c Hdf5ManagerCheetahLCLS.cpp
	g++ $(BEFORE) -c Hdf5ManagerCheetahSacla.cpp
	g++ $(BEFORE) -c Hdf5ManagerProcessing.cpp
	g++ $(BEFORE) -c Hdf5Table.cpp
	g++ $(BEFORE) -c Image.cpp
	g++ $(BEFORE) -c IndexManager.cpp
	g++ $(BEFORE) -c IndexingSolution.cpp
	g++ $(BEFORE) -c InputFileParser.cpp
	g++ $(BEFORE) -c LoggableObject.cpp
	g++ $(BEFORE) -c Logger.cpp
	g++ $(BEFORE) -c Masker.cpp
	g++ $(BEFORE) -c Matrix.cpp
	g++ $(BEFORE) -c Miller.cpp
	g++ $(BEFORE) -c MtzManager.cpp
	g++ $(BEFORE) -c MtzManagerRefine.cpp
	g++ $(BEFORE) -c MtzMerger.cpp
	g++ $(BEFORE) -c MtzRefiner.cpp
	g++ $(BEFORE) -c NelderMead.cpp
	g++ $(BEFORE) -c PNGFile.cpp
	g++ $(BEFORE) -c ProfileFit.cpp
	g++ $(BEFORE) -c PythonExt.cpp
	g++ $(BEFORE) -c RefinementGridSearch.cpp
	g++ $(BEFORE) -c RefinementStepSearch.cpp
	g++ $(BEFORE) -c RefinementStrategy.cpp
	g++ $(BEFORE) -c Reflection.cpp
	g++ $(BEFORE) -c Shoebox.cpp
	g++ $(BEFORE) -c SolventMask.cpp
	g++ $(BEFORE) -c SpectrumBeam.cpp
	g++ $(BEFORE) -c Spot.cpp
	g++ $(BEFORE) -c SpotFinder.cpp
	g++ $(BEFORE) -c SpotFinderCorrelation.cpp
	g++ $(BEFORE) -c SpotFinderQuick.cpp
	g++ $(BEFORE) -c SpotVector.cpp
	g++ $(BEFORE) -c StatisticsManager.cpp
	g++ $(BEFORE) -c TextManager.cpp
	g++ $(BEFORE) -c UnitCellLattice.cpp
	g++ $(BEFORE) -c Vector.cpp
	g++ $(BEFORE) -c main.cpp
	g++ $(BEFORE) -c misc.cpp
	
	g++ -g -o ../cppxfel *.o $(AFTER)
	
