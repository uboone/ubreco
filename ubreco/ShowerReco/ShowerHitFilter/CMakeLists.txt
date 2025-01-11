
art_make(	
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			#art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			#cetlib cetlib_except
                   	#larpandora_LArPandoraInterface	
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
			larcorealg::Geometry
			lardata::Utilities
			lardataobj::RecoBase
			larcore::Geometry_Geometry_service
          MODULE_LIBRARIES 
			ubreco::ShowerReco_ProximityClustering_Algorithms
		   	#larpandora_LArPandoraInterface	
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardataobj::RecoBase
			lardataobj::AnalysisBase	
			lardata::RecoObjects
                        #larreco::Deprecated  
			larreco::RecoAlg
			lardata::Utilities
                        ${Boost_SYSTEM_LIBRARY}
			larevt::Filters
			nusimdata::SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art::Persistency_Common 
                        #canvas
			art::Persistency_Provenance
			art::Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			#cetlib 
                        #cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
                        art_root_io::TFileService_service
        )

install_fhicl()
install_headers()
install_source()

