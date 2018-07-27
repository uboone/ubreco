art_make(
          LIB_LIBRARIES
			ubreco_ShowerReco_ClusterMerging_CMToolBase
			art_Framework_Core
			art_Framework_Services_Registry
			art_Framework_Principal
			lardataobj_RecoBase
			${ART_FRAMEWORK_CORE}
			larcorealg_Geometry	
			lardata_Utilities # added now
		   	#uboone_Geometry		
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

