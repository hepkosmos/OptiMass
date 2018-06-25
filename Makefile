####################################################################
# OptiMass Core SRC (libOptiMass), Header and Library information :
# libOptiMass.so/dylib should be installed in advance,
#     1) go to 'alm_base' directory
#     2) ./configure  (=> Check ROOT installation and MINUIT2 activation. Import env. variables)
#     3) make & make install (=> libOptiMass at 'lib'. headers at 'include/alm_base/')
LIBDIR   := lib
CXXFLAGS := -g -march=native -O2 -Wall -Wextra -std=c++11 -pedantic $(CXXFLAGS)
LDFLAGS  :=
LIBS     :=
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f

# alm_base
OPTIMASS := $(PWD)
CXXFLAGS += -I$(OPTIMASS)/include/alm_base
LIBS     += -Wl,-rpath,$(OPTIMASS)/lib -L$(OPTIMASS)/lib -lOptiMass

# ROOT
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --ldflags)
LIBS     += $(shell root-config --libs)
LIBS     += -lMinuit2

###################################################################################
# (*) ExRootAnalysis Information (also see 'README') :
# If your main.cpp does not depend on the ExRootAnalysis, just ignore this,
# however in order to use ExRootAnalysis as in the default skeleton main.cpp, then
# provide the path of your ExRootAnalysis directory to the varaible 'EXROOT' below.
EXROOT := $(EXROOT)
ifneq (${EXROOT},)
$(info )
$(info <= Your ExRootAnalysis directory = $(EXROOT) =>)
$(info )
else
$(info )
$(info <= Warning : ExRootAnalysis installation directory is not specified. =>)
$(info )
$(info     1> If your 'main.cpp' depends on ExRootAnalysis,)
$(info           provide your ExRootAnaysis installation directory to the varaible - 'EXROOT' in Makefile, otherwise followed compiling process will encounter errors from it. As for the installation of ExRootAnalysis, see 'README'.)
$(info )
$(info     2> If not, just ignore this message.)
$(info )
endif

# Model Card File Dictionary Information
MODELDIR := model/dict_src
MODELSRC := $(wildcard $(MODELDIR)/*.cpp)
MODELSRCh:= $(wildcard $(MODELDIR)/*.h)
MODELOBJ := $(MODELSRC:%.cpp=%.o)
MODELLIB := lib/libOptiMass-models.a
CXXFLAGS += -I$(MODELDIR)

# Main interface model
EXE    := optimass.x
EXESRC := main.cpp
EXEOBJ := $(EXESRC:%.cpp=%.o)

.PHONY: all clean model

all: $(EXE)
	@if test -f $(EXE);\
	then \
	echo "=======================================";\
	echo "Compilation finished. Run ./optimass.x";\
	echo "=======================================";\
	else \
	echo "========================================================";\
	echo "Compilation failed. Check your environment using README";\
	echo "========================================================";\
	fi

model: $(MODELLIB)

$(MODELLIB): CXXFLAGS += -fPIC
$(MODELLIB): $(MODELOBJ)
	$(MKDIR) $(LIBDIR)
	$(AR) $@ $^
	ranlib $@

$(EXE): CXXFLAGS += -I$(PWD) -I$(EXROOT)/ExRootAnalysis
$(EXE): LIBS     += -Wl,-rpath,$(EXROOT) -L$(EXROOT) -lExRootAnalysis
$(EXE): $(EXEOBJ) $(MODELOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) $(MODELSRC) $(MODELSRCh) $(MODELOBJ) $(MODELLIB)
	$(RM) $(EXEOBJ) $(EXE)

