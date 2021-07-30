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

		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}

        )

install_fhicl()
install_headers()
install_source()

