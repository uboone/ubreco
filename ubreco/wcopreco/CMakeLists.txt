
art_make(
	  LIB_LIBRARIES

			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}

			cetlib cetlib_except

          MODULE_LIBRARIES
			ubreco_wcopreco_data
			ubreco_wcopreco_deconv
			ubreco_wcopreco_algo
			ubreco_wcopreco_app
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}

        )

add_subdirectory(data)
add_subdirectory(deconv)
add_subdirectory(algo)
add_subdirectory(app)
add_subdirectory(dev)

install_headers()
#install_fhicl()
install_source()

