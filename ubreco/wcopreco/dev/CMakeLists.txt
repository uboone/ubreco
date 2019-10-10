include_directories(../data)
include_directories(../deconv)
include_directories(../algo)
include_directories(../app)

cet_make_exec(
	  main SOURCE main.cxx LIBRARIES
			${ROOT_BASIC_LIB_LIST}
			ubreco_wcopreco_data
			ubreco_wcopreco_deconv
			ubreco_wcopreco_algo
			ubreco_wcopreco_app
        )


#install_headers()
#install_fhicl()
#install_source()

