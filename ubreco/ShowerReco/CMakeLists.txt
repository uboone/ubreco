
art_make(
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
                   	#larpandora_LArPandoraInterface

          MODULE_LIBRARIES
		   	#larpandora_LArPandoraInterface
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
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
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
			CLHEP::CLHEP
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}

        )

add_subdirectory(TwoDimTools)
add_subdirectory(ShowerReco3D)
add_subdirectory(ClusterMerging)
add_subdirectory(ProximityClustering)
add_subdirectory(ReCluster)
add_subdirectory(Pi0Ana)
add_subdirectory(TruthFilters)

install_headers()
install_source()

