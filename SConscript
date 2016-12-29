import os
import libtbx.load_env
Import("env_etc")

env_etc.cppxfel_dist = libtbx.env.dist_path("cppxfel")
env_etc.cppxfel_include = os.path.dirname(env_etc.cppxfel_dist)
if (not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include")):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)
    env_etc.include_registry.append(
        env=env,
        paths=[
            env_etc.libtbx_include,
            env_etc.scitbx_include,
            env_etc.cctbx_include,
            env_etc.cctbx_include,
            env_etc.ccp4io_include,
            env_etc.boost_include,
            env_etc.boost_adaptbx_include,
            env_etc.python_include,
            env_etc.dxtbx_include,
            env_etc.cppxfel_include])
    print env_etc.cppxfel_dist
    env.Append(
		LIBS=env_etc.libm + ["scitbx_boost_python", "boost_thread", "boost_system",
		"boost_python",
		"cctbx",
		"ccp4io"])
    env.Replace(SHCCFLAGS=['-std=c++0x', '-fPIC', '-O3', '-g'])
    
if env_etc.clang_version:
  wd = ["-Wno-unused-variable"]
  env.Append(CCFLAGS=wd)

base_lib = libtbx.env.under_build(path="../base/lib")
base_include = libtbx.env.under_build(path="../base/include")
if (os.path.isdir(base_lib)):
  env.Append(LIBPATH=[base_lib])
  env.Append(CPPPATH=[base_include])

if 'BOOST_LOCATION' in os.environ:
	env.Append(LIBPATH = [os.environ['BOOST_LOCATION']])
	print ("Appending directory containing boost libraries: " + os.environ['BOOST_LOCATION'])

source = [
'source/Beam.cpp',
'source/GaussianBeam.cpp',
'source/SpectrumBeam.cpp',
'boost_python/cppxfel_ext.cc',
'source/AmbiguityBreaker.cpp',
'source/CSV.cpp',
'source/Detector.cpp',
'source/FileParser.cpp',
'source/FileReader.cpp',
'source/FreeLattice.cpp',
'source/FreeMillerLibrary.cpp',
'source/GeometryParser.cpp',
'source/GeometryRefiner.cpp',
'source/GraphDrawer.cpp',
'source/Hdf5Crystal.cpp',
'source/Hdf5Image.cpp',
'source/Hdf5Manager.cpp',
'source/Hdf5ManagerCheetah.cpp',
'source/Hdf5ManagerCheetahLCLS.cpp',
'source/Hdf5ManagerCheetahSacla.cpp',
'source/Hdf5ManagerProcessing.cpp',
'source/Hdf5Table.cpp',
'source/Image.cpp',
'source/IOMRefiner.cpp',
'source/IndexManager.cpp',
'source/IndexingSolution.cpp',
'source/InputFileParser.cpp',
'source/Logger.cpp',
'source/LoggableObject.cpp',
'source/Matrix.cpp',
'source/Miller.cpp',
'source/MtzGrouper.cpp',
'source/MtzMerger.cpp',
'source/MtzManager.cpp',
'source/MtzManagerMinimize.cpp',
'source/MtzManagerRefine.cpp',
'source/MtzRefiner.cpp',
'source/NelderMead.cpp',
'source/Panel.cpp',
'source/PanelParser.cpp',
'source/PNGFile.cpp',
'source/PythonExt.cpp',
'source/Reflection.cpp',
'source/RefinementStrategy.cpp',
'source/RefinementStepSearch.cpp',
'source/Shoebox.cpp',
'source/Spot.cpp',
'source/SpotVector.cpp',
'source/SpotFinder.cpp',
'source/SpotFinderQuick.cpp',
'source/SpotFinderCorrelation.cpp',
'source/SolventMask.cpp',
'source/StatisticsManager.cpp',
'source/UnitCellLattice.cpp',
'source/Vector.cpp',
'source/Wiki.cpp',
'source/main.cpp',
'source/misc.cpp']

env.SharedLibrary(
    target='#/lib/cppxfel_ext', 
    source=source,
    LIBS=env["LIBS"] + ['hdf5'] + ['hdf5_hl'] + ['png'])
