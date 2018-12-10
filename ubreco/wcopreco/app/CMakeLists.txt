include_directories(../data)
include_directories(../deconv)
include_directories(../algo)
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(
	  LIB_LIBRARIES
			${ROOT_BASIC_LIB_LIST}
			ubreco_wcopreco_data
			ubreco_wcopreco_deconv
			ubreco_wcopreco_algo
        )


install_headers()
#install_fhicl()
install_source()

