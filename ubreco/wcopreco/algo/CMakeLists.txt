include_directories(../data)
include_directories(../deconv)
include_directories($ENV{EIGEN_INC})

art_make(
          LIB_LIBRARIES

                        ubreco_wcopreco_data
                        ubreco_wcopreco_deconv
        )


install_headers()
#install_fhicl()
install_source()
