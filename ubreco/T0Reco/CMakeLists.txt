cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH ) 
include_directories ( $ENV{CRYHOME}/src )

art_make(
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}

			cetlib cetlib_except

          MODULE_LIBRARIES
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			larcoreobj_SummaryData
			larevt_Filters
			nusimdata_SimulationBase
			nurandom_RandomUtils_NuRandomService_service
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
                        ubevt_Utilities						
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
			${ROOT_PHYSICS}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_CORE}
			${ROOT_EG}
			${CRY}
        )

install_headers()
install_fhicl()
install_source()

