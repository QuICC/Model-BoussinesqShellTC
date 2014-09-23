#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%% Create all the test targets %%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#
# Install MatrixMarket test files
#
INSTALL(DIRECTORY ../TestSuite/SparseLinearSolver DESTINATION
   ${PROJECT_BINARY_DIR}/TestSuite
   FILES_MATCHING PATTERN "*.mtx")

add_executable("SparseLinearSolverTest" "SparseLinearSolverTest.cpp"
   "../src/Framework/FrameworkBase.cpp"
   "../src/Framework/SerialFramework.cpp"
   "../src/Exceptions/Exception.cpp"
   "../src/Timers/ITimer.cpp"
   "../src/Timers/SerialTimer.cpp"
   "../src/Profiler/ProfilerBase.cpp"
   "../src/Profiler/SerialProfiler.cpp"
   )
