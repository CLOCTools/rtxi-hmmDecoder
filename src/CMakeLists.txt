set(LIB_SOURCES ${LIB_SOURCES})
file(READ sourcelist.txt LIB_SOURCES)
message(STATUS "LIB_SOURCES           = ${LIB_SOURCES}")

target_sources(${CMAKE_PROJECT_NAME} PRIVATE ${LIB_SOURCES})