SET(NAME Sources)
SET(SOURCES initialisation.cpp inlet.cpp outlet.cpp
setGhostCells.cpp setting.cpp timeStep.cpp wall.cpp)

ADD_LIBRARY(${NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${NAME} Fvm)
