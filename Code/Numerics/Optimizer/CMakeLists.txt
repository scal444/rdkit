
rdkit_library(Optimizer
              BFGSOpt.cpp BFGSKernels.cpp LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib Trajectory RDGeneral)
target_compile_definitions(Optimizer PRIVATE RDKIT_OPTIMIZER_BUILD)

rdkit_headers(BFGSOpt.h BFGSKernels.h DEST Numerics/Optimizer)

rdkit_catch_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer )


