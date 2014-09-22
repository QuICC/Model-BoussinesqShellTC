#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%% Create all the test targets %%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

add_executable("SparseLinearSolverTest" "SparseLinearSolverTest.cpp"
   "../src/Framework/FrameworkBase.cpp"
   "../src/Framework/SerialFramework.cpp"
   "../src/Exceptions/Exception.cpp"
   "../src/Profiler/ProfilerBase.cpp"
   "../src/Profiler/SerialProfiler.cpp"
   )
