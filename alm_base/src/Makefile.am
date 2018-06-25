AM_CXXFLAGS = -std=c++11 -I$(top_srcdir)/src $(ROOTCFLAGS)

if DEBUG
AM_CXXFLAGS += -O0 -Wall -Wextra -pedantic
endif

lib_LTLIBRARIES        = libOptiMass.la
libOptiMass_la_SOURCES = ALMController.cpp MassFunction.cpp \
	MassFunctionInterface.cpp MassMinimizer.cpp MathUtils.cpp \
	ProcessTree.cpp StringUtils.cpp
libOptiMass_la_LIBADD  = -L$(ROOTLIBDIR) $(ROOTLIBS) -lMinuit2

pkginclude_HEADERS = ALMController.h ConstraintBase.h MassFunction.h \
	MassFunctionInterface.h MassMinimizer.h MathUtils.h ProcessTree.h \
	StringUtils.h Types.h
