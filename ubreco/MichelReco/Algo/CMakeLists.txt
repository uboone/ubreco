art_make(
          LIB_LIBRARIES
                        ubreco_MichelReco_Fmwk
                        #messagefacility::MF_MessageLogger
                        #
                        #cetlib::cetlib
                        ROOT::Graf
                        ROOT::Spectrum
        )


install_headers()
install_fhicl()
install_source()
