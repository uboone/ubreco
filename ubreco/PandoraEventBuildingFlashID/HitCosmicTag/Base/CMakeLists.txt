set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} ${ROOTLIB} )
art_make( BASENAME_ONLY
   	      LIBRARY_NAME     ubreco_PandoraEventBuildingFlashID_HitCosmicTag_Base
	      LIB_LIBRARIES
	      ${ROOT_BASIC_LIB_LIST}
	      fhiclcpp::fhiclcpp
              cetlib::cetlib cetlib_except
	    )
install_headers()
install_source()
