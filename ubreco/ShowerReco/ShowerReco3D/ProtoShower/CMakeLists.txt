art_make(
          LIB_LIBRARIES
			art_Framework_Core
			art_Framework_Services_Registry
			art_Framework_Principal
			art_Persistency_Provenance
			${ART_FRAMEWORK_CORE}
			larcorealg_Geometry	
			lardataobj_RecoBase
		   	larcorealg_Geometry		
			larcore_Geometry_Geometry_service	
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB} cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}

        )


install_headers()
install_source()

