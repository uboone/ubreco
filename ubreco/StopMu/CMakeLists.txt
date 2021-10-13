include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
#cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH )
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories ( $ENV{CRYHOME}/src )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make(	
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
                   	#larpandora_LArPandoraInterface	
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
			larcorealg_Geometry
			lardata_Utilities
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
			lardataobj_MCBase
			larcore_Geometry_Geometry_service
          MODULE_LIBRARIES 
		   	#larpandora_LArPandoraInterface	
			larevt_CalibrationDBI_IOVData
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
			lardataobj_AnalysisBase	
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			nurandom_RandomUtils_NuRandomService_service
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
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
			CLHEP::CLHEP
			${ROOT_PHYSICS}
      			${ROOT_GEOM}
      			${ROOT_XMLIO}
      			${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_CORE}
			${ROOT_EG}
			CRY
        )


install_fhicl()
install_headers()
install_source()
