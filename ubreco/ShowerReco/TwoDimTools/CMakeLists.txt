
art_make(	
	  LIB_LIBRARIES
			art_Framework_Core
			art_Framework_Services_Registry
			art_Framework_Principal ${ART_FRAMEWORK_CORE}
			larcorealg_Geometry lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			art_Utilities canvas ${MF_MESSAGELOGGER}
			${MF_UTILITIES} ${CETLIB}
			${Boost_SYSTEM_LIBRARY} ${ROOT_BASIC_LIB_LIST}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${MF_MESSAGELOGGER}
			cetlib cetlib_except

        )

install_headers()
install_source()

