art_make( BASENAME_ONLY
    	  LIBRARY_NAME     ubreco_SelectionTool_OpT0FinderAlgorithms
	  LIB_LIBRARIES
                	  larsim_PhotonPropagation
			  larsim_PhotonPropagation_PhotonVisibilityService_service
			  larcorealg_Geometry
			  ubcore::ubana_BasicTool_GeoAlgo
			  ubreco_SelectionTool_OpT0FinderBase
			  ubobj_Trigger
                          UBOpReadoutMap_service
  		          art::Framework_Core
   		          art::Framework_Principal
                          art::Framework_Services_Registry
			  art::Persistency_Common
			  art::Persistency_Provenance
			  art::Utilities
                          canvas::canvas
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib
                          cetlib_except::cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
