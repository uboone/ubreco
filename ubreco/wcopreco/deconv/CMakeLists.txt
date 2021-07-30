include_directories(../data)
include_directories($ENV{EIGEN_INC})

art_make(
          LIB_LIBRARIES
                        ROOT::Hist
                        ubreco_wcopreco_data
        )


install_headers()
#install_fhicl()
install_source()
