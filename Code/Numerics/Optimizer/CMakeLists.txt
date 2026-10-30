
rdkit_select_kernels(BFGS_KERNELS BFGSKernels.cpp BFGSKernelsHwy.cpp)
rdkit_library(Optimizer
              BFGSOpt.cpp ${BFGS_KERNELS} LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib Trajectory RDGeneral)
target_compile_definitions(Optimizer PRIVATE RDKIT_OPTIMIZER_BUILD)
rdkit_link_highway(Optimizer)

rdkit_headers(BFGSOpt.h BFGSKernels.h DEST Numerics/Optimizer)

rdkit_catch_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer )
rdkit_link_highway(testOptimizer)
