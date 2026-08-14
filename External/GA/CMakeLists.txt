rdkit_library(ga ga/StringChromosome.cpp ga/Chromosome.cpp
                 ga/IntegerStringChromosomePolicy.cpp ga/BinaryStringChromosomePolicy.cpp util/RandomUtil.cpp util/Util.cpp)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND
   CMAKE_CXX_FLAGS MATCHES "local-bounds")
  target_compile_options(ga PRIVATE
    "-fsanitize-ignorelist=${CMAKE_SOURCE_DIR}/Code/Sanitizers/ubsan.ignorelist")
endif()
