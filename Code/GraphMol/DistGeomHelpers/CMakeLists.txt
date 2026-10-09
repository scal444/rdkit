
rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp EmbedderUtils.cpp ZMatrixBuilder.cpp
              LINK_LIBRARIES MolAlign ForceFieldHelpers SubstructMatch GraphMol DistGeometry Alignment
                )
target_compile_definitions(DistGeomHelpers PRIVATE RDKIT_DISTGEOMHELPERS_BUILD)

rdkit_headers(BoundsMatrixBuilder.h BoundsMatrixBuilderDetails.h ZMatrixBuilder.h
              Embedder.h DEST GraphMol/DistGeomHelpers)

rdkit_catch_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers MolTransforms FileParsers SmilesParse SHARD_COUNT 4 )

rdkit_catch_test(distGeomHelpersCatch catch_tests.cpp 
LINK_LIBRARIES DistGeomHelpers MolTransforms FileParsers SmilesParse SHARD_COUNT 4 )

rdkit_catch_test(testDgeomHelpersGH3019 testGithub3019.cpp
           LINK_LIBRARIES DistGeomHelpers SmilesParse SHARD_COUNT 2)


if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()

# Compile public consumers without any source-tree RDKit include directories.
if(RDK_BUILD_CPP_TESTS AND BUILD_TESTING AND RDK_INSTALL_DEV_COMPONENT
   AND NOT RDK_INSTALL_INTREE)
  set(installed_header_test_prefix "${CMAKE_CURRENT_BINARY_DIR}/installed-header-test")
  set(installed_header_test_sources)
  set(installed_header_test_inputs
      testInstalledBoundsMatrixBuilder.cpp testInstalledTorsionPreferences.cpp)
  if(RDK_BUILD_THREADSAFE_SSS)
    list(APPEND installed_header_test_inputs
         testInstalledSynthonSpaceHitSet.cpp testInstalledSynthonSpaceShapeSearcher.cpp)
  endif()
  foreach(source IN LISTS installed_header_test_inputs)
    configure_file(${source} ${CMAKE_CURRENT_BINARY_DIR}/${source} COPYONLY)
    list(APPEND installed_header_test_sources ${CMAKE_CURRENT_BINARY_DIR}/${source})
  endforeach()
  add_library(testInstalledHeaders OBJECT EXCLUDE_FROM_ALL
              ${installed_header_test_sources})
  # Replace inherited include directories; rdkit_base would reintroduce Code/.
  set_property(TARGET testInstalledHeaders PROPERTY INCLUDE_DIRECTORIES
               "${installed_header_test_prefix}/${RDKit_HdrDir}")
  target_link_libraries(testInstalledHeaders PRIVATE Boost::boost)
  if(TARGET Eigen3::Eigen)
    target_link_libraries(testInstalledHeaders PRIVATE Eigen3::Eigen)
  endif()
  configure_file(testInstalledHeaders.cmake.in testInstalledHeaders.cmake @ONLY)
  add_test(NAME testInstalledHeaders
           COMMAND ${CMAKE_COMMAND} -DTEST_CONFIG=$<CONFIG>
                   -P ${CMAKE_CURRENT_BINARY_DIR}/testInstalledHeaders.cmake)
  set_tests_properties(testInstalledHeaders PROPERTIES RUN_SERIAL TRUE)
endif()
