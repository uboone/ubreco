set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} ${ROOTLIB} )
art_make( BASENAME_ONLY
   	  LIBRARY_NAME     ubreco_SelectionTool_OpT0FinderBase
	  LIB_LIBRARIES
	  ubana_BasicTool_GeoAlgo
	  ${ROOT_BASIC_LIB_LIST}
	  fhiclcpp::fhiclcpp
	  art::Utilities
          cetlib::cetlib cetlib_except

	larcorealg_Geometry
	larcore_Geometry_Geometry_service
	lardataalg_DetectorInfo
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Persistency_Provenance
	art::Framework_Services_Registry
	canvas::canvas
	)
install_headers()
install_source()
