include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

add_subdirectory(job)
add_subdirectory(HitCosmicTag)


art_make(         TOOL_LIBRARIES
		  	larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardataobj_AnalysisBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        larreco_RecoAlg_ClusterRecoUtil
                        larreco_RecoAlg
                        larreco_ClusterFinder
                        larreco_Calorimetry
                        ubreco_PandoraEventBuildingFlashID_HitCosmicTag_Base
                        ubreco_PandoraEventBuildingFlashID_HitCosmicTag_Algorithms
			PandoraSDK
                        PandoraMonitoring
                        larpandora_LArPandoraInterface
                        larpandora_LArPandoraEventBuilding
                        nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        
                        Boost::filesystem
                        ROOT::Geom
                        ubreco_SelectionTool_OpT0FinderBase
                        ubreco_SelectionTool_OpT0FinderAlgorithms
                        ubevt_Utilities

          	MODULE_LIBRARIES
	  		ubreco_ShowerReco_ShowerReco3D_Base
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_AnalysisBase
			lardataobj_Simulation
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
                        larpandora_LArPandoraInterface
                        larpandora_LArPandoraEventBuilding
			larevt_Filters
			nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Framework_Services_Optional_RandomNumberGenerator_service
			
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ubreco_SelectionTool_OpT0FinderBase
                        ubreco_SelectionTool_OpT0FinderAlgorithms
                        ubevt_Utilities
)

install_headers()
install_fhicl()
install_source()
