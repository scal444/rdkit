add_executable(bench EXCLUDE_FROM_ALL bench_common.cpp
  descriptors.cpp
  fingerprint.cpp
  meta.cpp
  mol.cpp
  molops.cpp
  pickle.cpp
  smiles.cpp
  stereo.cpp
  substruct_match.cpp
)
target_link_libraries(bench rdkitCatch
  CIPLabeler
  Descriptors
  Fingerprints
  SmilesParse
)

add_executable(smiles_pipeline_bench EXCLUDE_FROM_ALL smiles_pipeline_bench.cpp)
target_link_libraries(smiles_pipeline_bench SmilesParse)

add_executable(chembl_ops_bench EXCLUDE_FROM_ALL chembl_ops_bench.cpp)
target_link_libraries(chembl_ops_bench
  DataStructs
  Descriptors
  DistGeomHelpers
  Fingerprints
  GraphMol
  SmilesParse
)
if(RDK_BUILD_INCHI_SUPPORT)
  target_compile_definitions(chembl_ops_bench PRIVATE RDK_BUILD_INCHI_SUPPORT)
  target_include_directories(chembl_ops_bench PRIVATE
    ${CMAKE_SOURCE_DIR}/External
  )
  target_link_libraries(chembl_ops_bench RDInchiLib)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(chembl_ops_bench PRIVATE
    -fno-lto
    -fno-profile-instr-generate
    -fno-profile-instr-use
  )
  target_link_options(chembl_ops_bench PRIVATE
    -fno-lto
    -fno-profile-instr-generate
    -fno-profile-instr-use
  )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(chembl_ops_bench PRIVATE
    -fno-lto
    -fno-profile-generate
    -fno-profile-use
  )
  target_link_options(chembl_ops_bench PRIVATE
    -fno-lto
    -fno-profile-generate
    -fno-profile-use
  )
endif()

if(RDK_BUILD_INCHI_SUPPORT)
  target_link_libraries(bench RDInchiLib)
  target_sources(bench PRIVATE inchi.cpp)
endif(RDK_BUILD_INCHI_SUPPORT)

if(RDK_BUILD_CPP_TESTS)
  # add a fast version of the benchmarks to the default unit tests
  # to protect the benchmarks from bit-rot
  add_test(
    NAME quickbench
    COMMAND bench --benchmark-samples 1 --benchmark-warmup-time 0
  )
  
  # ctest does not automatically build test executables
  # so, like with other unit tests in this repo, we need to manually
  # build the target before running ctest. This usually gets done with
  # an `all` or an `install` build. We re-add bench to he `all` group
  # to keep it similar to the other unit tests
  set_target_properties(bench PROPERTIES EXCLUDE_FROM_ALL FALSE)
endif(RDK_BUILD_CPP_TESTS)
