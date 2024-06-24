file(
	GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(
	Test_Decompositions
	${SRC_FILES}
)

#target_link_libraries(
#	Test_Decompositions PRIVATE
#	Decompositions
#)

target_include_directories(
	Test_Decompositions PUBLIC 
	${CMAKE_SOURCE_DIR}/include
)

target_compile_features(Test_Decompositions PRIVATE cxx_std_17)