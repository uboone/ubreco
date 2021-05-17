
art_make(	
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			nutools_ParticleNavigation
			lardataobj_Simulation
			lardataobj_RawData
			larsim_Simulation
			lardataobj_AnalysisBase	
			lardata_RecoObjects 
			larreco_Deprecated  
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			SignalShapingServiceMicroBooNE_service
			ubobj_Trigger
			ubobj_Optical
			${LARRECO_LIB}
			${LARDATA_LIB}
            ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${FHICLCPP}
			cetlib cetlib_except
			cetlib 
			${CLHEP}
            ${ROOT_GEOM}
            ${ROOT_XMLIO}
            ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )

add_subdirectory(Analysis)
add_subdirectory(HitFinding)
add_subdirectory(CosmicBackground)

install_headers()
install_fhicl()
install_source()

