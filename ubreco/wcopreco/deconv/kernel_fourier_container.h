#ifndef KERNEL_FOURIER_CONTAINER_H
#define KERNEL_FOURIER_CONTAINER_H



#include "kernel_fourier.h"
#include <vector>

namespace wcopreco {

  class kernel_fourier_container : public std::vector<kernel_fourier*> {

  public:
    kernel_fourier_container();
    virtual ~kernel_fourier_container() ;


    void add_kernel(kernel_fourier *kernel) { push_back(kernel);}
    void clear_kernels_v() {clear();}


  protected:

  };

}

#endif
