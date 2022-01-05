
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          MODULE_LIBRARIES  larsim_Simulation 
			    lardataobj_Simulation
			    lardata_RecoObjects
			    lardataobj_RecoBase
			    larcore_Geometry_Geometry_service
			    larcorealg_GeoAlgo
			    larcorealg_GeoAlgo_dict
			    larcorealg_Geometry
			    lardataobj_AnalysisBase
			    lardata_Utilities
			    nusimdata::SimulationBase
			    art::Framework_Core
			    art::Framework_Principal
			    art::Framework_Services_Registry
			    
			    art::Framework_Services_Optional_RandomNumberGenerator_service
			    art::Persistency_Common canvas
			    art::Persistency_Provenance canvas
			    art::Utilities
			    messagefacility::MF_MessageLogger
			    
			    fhiclcpp::fhiclcpp
			    CLHEP::CLHEP
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
			    larreco_RecoAlg
			    larreco_RecoAlg_ClusterRecoUtil
			    larreco_ClusterFinder
			    ${ROOT_GENVECTOR}
         )

install_headers()
install_fhicl()
install_source()


