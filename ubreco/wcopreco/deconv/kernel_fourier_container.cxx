#include "kernel_fourier_container.h"
#include <vector>

using namespace wcopreco ;

wcopreco::kernel_fourier_container::kernel_fourier_container()
  {

  }

wcopreco::kernel_fourier_container::~kernel_fourier_container()
  {
    for (size_t i=0; i<size(); i++){
      delete at(i);
    }
  }
