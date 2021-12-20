
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
#cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH )
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories ( $ENV{CRYHOME}/src )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make( BASENAME_ONLY
      MODULE_LIBRARIES
			larevt_CalibrationDBI_IOVData
      			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			lardataobj_RecoBase
			lardata_RecoObjects
			larreco_RecoAlg
			ubobj_MuCS
			ubevt::SignalShapingServiceMicroBooNE_service
			lardataobj_AnalysisBase
      			nusimdata::SimulationBase
      			${LARRECO_LIB}
      			${LARDATA_LIB}
      			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art::Persistency_Common canvas::canvas
			art::Persistency_Provenance canvas::canvas
			art::Utilities canvas::canvas
			nurandom::RandomUtils_NuRandomService_service
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
			CLHEP::CLHEP
			ROOT::Physics
      			ROOT::Geom
      			ROOT::XMLIO
      			ROOT::Gdml

			ROOT::Core
			ROOT::EG
			SQLite::SQLite3
			CRY
			IFDH
)

install_fhicl()
install_source()
install_headers()
