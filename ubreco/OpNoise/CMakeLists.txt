
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
			    nug4_MagneticField_MagneticField_service
			    nusimdata_SimulationBase
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    art_Persistency_Common canvas
			    art_Persistency_Provenance canvas
			    art_Utilities
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
			    ${FHICLCPP}
			    ${CLHEP}
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


