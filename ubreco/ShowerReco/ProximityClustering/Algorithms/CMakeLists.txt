art_make(
          LIB_LIBRARIES
		   	larcorealg_Geometry		
			lardata_Utilities	
			larcore_Geometry_Geometry_service	
			lardataobj_RecoBase
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB} cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}

        )

install_fhicl()
install_headers()
install_source()

