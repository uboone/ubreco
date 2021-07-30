art_make(
          LIB_LIBRARIES
                        ubreco_MichelReco_Fmwk
                        #${MF_MESSAGELOGGER}
                        #
                        #${CETLIB}
                        ROOT::Graf
                        ${ROOT_SPECTRUM}
        )


install_headers()
install_fhicl()
install_source()
