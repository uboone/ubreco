include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

cet_make( LIBRARIES ${ROOT_BASIC_LIB_LIST}
	            ${ROOT_EG}
	            NO_DICTIONARY
	)

art_dictionary( DICTIONARY_LIBRARIES ubreco_WcpPortedReco_ProducePort_WCP )

install_headers()
install_fhicl()
install_source()
