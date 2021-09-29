art_make(
          LIB_LIBRARIES
		   	larcorealg_Geometry
			lardata_Utilities
			larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
			lardataobj_RecoBase
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB} cetlib_except

		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}

        )

install_fhicl()
install_headers()
install_source()

