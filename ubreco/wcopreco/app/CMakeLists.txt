include_directories( $ENV{EIGEN_INC} )
include_directories(../data)
include_directories(../deconv)
include_directories(../algo)
include_directories($ENV{EIGEN_INC})

art_make(
          LIB_LIBRARIES

                        ubreco_wcopreco_data
                        ubreco_wcopreco_deconv
                        ubreco_wcopreco_algo
                        ROOT::Core
        )


install_headers()
#install_fhicl()
install_source()
