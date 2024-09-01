file(
	GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../include/*.hpp
)

add_executable(
	Test_Decompositions
	${SRC_FILES}
)

target_include_directories(
	Test_Decompositions PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_compile_features(Test_Decompositions PRIVATE cxx_std_17)