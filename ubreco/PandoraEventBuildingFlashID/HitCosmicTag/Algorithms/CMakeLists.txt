# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} ${ROOTLIB} )

art_make( BASENAME_ONLY
    	    LIBRARY_NAME     ubreco_PandoraEventBuildingFlashID_HitCosmicTag_Algorithms
	        LIB_LIBRARIES
                           ubreco_PandoraEventBuildingFlashID_HitCosmicTag_Base
			   art::Framework_Core
   		           art::Framework_Principal
                           art::Framework_Services_Registry
     		           
     		           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           ${ROOT_BASIC_LIB_LIST}
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib cetlib_except
	      )

install_headers()
install_source()
