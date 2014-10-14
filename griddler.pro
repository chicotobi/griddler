TEMPLATE = app
CONFIG  -= qt
CONFIG += debug
LIBS	+= -lm
QMAKE_CXXFLAGS   += -std=c++11 -fopenmp
QMAKE_CFLAGS   += -std=c++11 -fopenmp
QMAKE_LFLAGS   += -std=c++11 -fopenmp
SOURCES  = main.cpp
TARGET   = Griddler
