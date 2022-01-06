art_make(

MODULE_LIBRARIES 
	ubobj_WcpPort
	larcorealg_Geometry
	larcore_Geometry_Geometry_service
	lardataobj_RawData
	lardataobj_RecoBase
	lardataobj_MCBase
	nusimdata::SimulationBase
	lardata_Utilities
	lardataalg_DetectorInfo
	larreco_RecoAlg
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Persistency_Provenance
	art::Utilities
	canvas::canvas
	messagefacility::MF_MessageLogger
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	ROOT::EG
)

#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory(Check)
add_subdirectory(job)

install_fhicl()
install_headers()
install_source()
