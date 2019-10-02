art_make( MODULE_LIBRARIES larcorealg_Geometry
	larcore_Geometry_Geometry_service
	lardataobj_RawData
	lardataobj_RecoBase
	lardataobj_MCBase
	nusimdata_SimulationBase
	lardata_Utilities
	lardataalg_DetectorInfo
	larreco_RecoAlg
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_PRINCIPAL}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_FRAMEWORK_SERVICES_OPTIONAL}
	${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	art_Persistency_Common
	art_Persistency_Provenance
	art_Utilities
	canvas
	${MF_MESSAGELOGGER}

	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_GEOM}
	${ROOT_XMLIO}
	${ROOT_GDML}
	${ROOT_EG}
	${ROOT_BASIC_LIB_LIST}
)

install_fhicl()
install_headers()
install_source()
