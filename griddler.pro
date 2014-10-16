TEMPLATE = app
CONFIG  -= qt
CONFIG += debug
debug {
QMAKE_CXXFLAGS   += -std=c++11 -pg
QMAKE_CFLAGS   += -std=c++11 -pg
QMAKE_LFLAGS   += -std=c++11 -pg
}
!debug {
QMAKE_CXXFLAGS   += -std=c++11 -O2 -fopenmp
QMAKE_CFLAGS   += -std=c++11 -O2 -fopenmp
QMAKE_LFLAGS   += -std=c++11 -O2 -fopenmp
}

SOURCES  = main.cpp
TARGET   = Griddler
INCLUDEPATH += /u/h/hofmannt/pHofmannT/boost/boost_1_56_0
