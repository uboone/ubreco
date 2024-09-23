CPPFLAGS=-I $(BOOST_INC) \
	 -I $(CANVAS_INC) \
	 -I $(CANVAS_ROOT_IO_INC) \
	 -I $(CETLIB_INC) \
	 -I $(CETLIB_EXCEPT_INC) \
	 -I $(FHICLCPP_INC) \
	 -I $(GALLERY_INC) \
	 -I $(LARCOREOBJ_INC) \
	 -I $(LARDATAOBJ_INC) \
	 -I $(LARCOREALG_INC) \
	 -I $(LARPANDORA_INC) \
	 -I $(LARDATA_INC) \
	 -I $(LARRECO_INC) \
	 -I $(LAREVT_INC) \
	 -I $(NUSIMDATA_INC) \
	 -I $(ROOT_INC) \
	 -I $(LARSIM_INC) \
	 -I $(CLHEP_INC) \
	 -I $(GENIE_INC) \
	 -I $(LARSIM_INC) \
	 -I $(HEP_CONCURRENCY_INC)

CXX=g++
CXXFLAGS=-std=c++17 -Wall -pedantic -g -I.

LDFLAGS=$(shell root-config --libs) \
	-L $(CANVAS_LIB) -l canvas \
	-L $(CETLIB_LIB) -l cetlib \
        -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
	-L $(GALLERY_LIB) -l gallery \
	-L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
	-L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
	-L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase \
	-L $(LARSIM_LIB) -l larsim_EventWeight_Base  \
	-L $(LARRECO_LIB) -l larreco_RecoAlg \
	-L $(FHICLCPP_LIB) -l fhiclcpp

all:  clippy

dumpwgh: clippy.cpp
	@echo Building $@
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

