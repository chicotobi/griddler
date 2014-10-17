TEMPLATE = app
CONFIG  -= qt
CONFIG -= debug
CONFIG += openmp
debug {
QMAKE_CXXFLAGS   += -std=c++11 -pg
QMAKE_CFLAGS   += -std=c++11 -pg
QMAKE_LFLAGS   += -std=c++11 -pg
}
!debug {
QMAKE_CXXFLAGS   += -std=c++11 -O2
QMAKE_CFLAGS   += -std=c++11 -O2
QMAKE_LFLAGS   += -std=c++11 -O2 
}
openmp {
QMAKE_CXXFLAGS   += -fopenmp
QMAKE_CFLAGS   += -fopenmp
QMAKE_LFLAGS   += -fopenmp
}

SOURCES  = main.cpp
TARGET   = Griddler
INCLUDEPATH += /u/h/hofmannt/pHofmannT/boost/boost_1_56_0
