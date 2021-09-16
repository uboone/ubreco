include_directories( $ENV{EIGEN_INC} )
include_directories(../data)
include_directories(../deconv)

art_make(
	  LIB_LIBRARIES
			${ROOT_BASIC_LIB_LIST}
			ubreco_wcopreco_data
			ubreco_wcopreco_deconv
        )


install_headers()
#install_fhicl()
install_source()

