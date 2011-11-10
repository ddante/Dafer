cmake_minimum_required (VERSION 2.6)
project (RD_Sca Fortran)

message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")

SET(EXECUTABLE_OUTPUT_PATH ../build)

SET(SOURCES Src/Elements/Element_Class.f90
            Src/Elements/Face_Class.f90
            Src/Elements/Segment_Class.f90
            Src/Elements/Triangle_Class.f90
	    Src/Elements/Quadrangle_Class.f90
	    Src/Init_problem.f90
            Src/Models.f90
	    Src/Geometry.f90
	    Src/Time_integration.f90
	    Src/Space_integration.f90
            Src/Main.f90
	    Src/Post_pro.f90
   )

ADD_EXECUTABLE(solver.exe ${SOURCES})


