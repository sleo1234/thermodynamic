project(thermo)
set(CMAKE_CXX_STANDARD 17)
SET(SOURCE_FILES Solver.cpp MatrixOperations.cpp PropertyPackage.cpp ExpressionBuilder.cpp FlashCalculation.cpp main.cpp)
add_executable(Solver.cpp MatrixOperations.cpp PropertyPackage.cpp ExpressionBuilder.cpp FlashCalculation.cpp main.cpp)