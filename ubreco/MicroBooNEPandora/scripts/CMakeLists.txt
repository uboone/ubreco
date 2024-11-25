
# The xml files here are used to initialize pandora
# They need to be found via FW_SEARCH_PATH

# install xml files
file(GLOB xml_file *.xml)
file(GLOB pt_file *.pt)
install_fw(LIST ${xml_file})

if(pt_file)
    install_fw(LIST ${pt_file})
endif()
