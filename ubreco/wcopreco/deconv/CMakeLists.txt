include_directories(../data)

art_make(
	  LIB_LIBRARIES
			${ROOT_BASIC_LIB_LIST}
			ubreco_wcopreco_data
        )


install_headers()
#install_fhicl()
install_source()

