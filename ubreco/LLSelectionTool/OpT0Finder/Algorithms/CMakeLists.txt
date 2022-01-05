# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
    	  LIBRARY_NAME     ubreco_SelectionTool_OpT0FinderAlgorithms
	  LIB_LIBRARIES
                	  larsim_PhotonPropagation
			  larsim_PhotonPropagation_PhotonVisibilityService_service
			  larcorealg_Geometry
			  ubana_BasicTool_GeoAlgo
			  ubreco_SelectionTool_OpT0FinderBase
			  ubobj_Trigger
                          UBOpReadoutMap_service
  		          art::Framework_Core
   		          art::Framework_Principal
                          art::Framework_Services_Registry
     		          
			  art::Persistency_Common canvas
			  art::Persistency_Provenance canvas
			  art::Utilities canvas
                          ${ROOT_BASIC_LIB_LIST}
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
