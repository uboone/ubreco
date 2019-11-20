
# where should the scripts/..xml file be installed?  Perhaps in bin?

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

set(TORCH_LIBRARY1 "$ENV{LIBTORCH_LIB}/libtorch.so")
set(TORCH_LIBRARY2 "$ENV{LIBTORCH_LIB}/libcaffe2.so")
set(TORCH_LIBRARY3 "$ENV{LIBTORCH_LIB}/libc10.so")
set(TORCH_INCLUDE_DIR "$ENV{LIBTORCH_INC}")
link_libraries(${TORCH_LIBRARY1} ${TORCH_LIBRARY2} ${TORCH_LIBRARY3})
include_directories(${TORCH_INCLUDE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  MicroBooNEPandora
          LIB_LIBRARIES larpandora_LArPandoraInterface
                        ${LARPANDORACONTENT}
                        ${PANDORASDK}
                        ${PANDORAMONITORING}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
                        ${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES MicroBooNEPandora
          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)

