SET(NAME Saving)
SET(SOURCES saveResults.cpp saveNormRez.cpp vtk_output.cpp)

ADD_LIBRARY(${NAME} ${SOURCES})
