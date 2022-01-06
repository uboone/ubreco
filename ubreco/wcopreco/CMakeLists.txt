
art_make(
	  LIB_LIBRARIES

			art::Framework_Core
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger

			cetlib::cetlib
                        cetlib_except::cetlib_except

          MODULE_LIBRARIES
			ubreco_wcopreco_data
			ubreco_wcopreco_deconv
			ubreco_wcopreco_algo
			ubreco_wcopreco_app
			cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml

        )

add_subdirectory(data)
add_subdirectory(deconv)
add_subdirectory(algo)
add_subdirectory(app)
add_subdirectory(dev)

install_headers()
#install_fhicl()
install_source()

