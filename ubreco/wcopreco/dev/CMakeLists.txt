include_directories(../data)
include_directories(../deconv)
include_directories(../algo)
include_directories(../app)

cet_make_exec(
  NAME wcopreco
  SOURCE main.cxx
  LIBRARIES
  PUBLIC
  ubreco::wcopreco_app
  Eigen3::Eigen
  ROOT::Tree
)

install_headers()
#install_fhicl()
install_source()
