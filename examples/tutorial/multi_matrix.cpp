/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#define VIENNACL_PROFILING_ENABLE
//#define VIENNACL_DEBUG_SCHEDULER
//#define VIENNACL_DEBUG_ALL

//#define VIENNACL_WITH_OPENCL

// include necessary system headers
#include <iostream>

#include "CL/cl.h"

#include "viennacl/distributed/multi_matrix.hpp"
//#include "viennacl/distributed/multi_vector.hpp"

#include "viennacl/distributed/scheduler.hpp"
#include "viennacl/distributed/utils.hpp"

#include "viennacl/linalg/prod.hpp"

#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/event.hpp"


#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

#define N_SIZES 22
#define SIZE_INC 1024



//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

template<typename ScalarType>
double run_benchmark(unsigned int size, double & vm_usage, double & rss){


    typedef viennacl::distributed::multi_matrix<ScalarType> gpu_mat_t;


    gpu_mat_t C(size,size);
    gpu_mat_t A(size,size);
    gpu_mat_t B(size,size);

    viennacl::distributed::timer t;

//    gpu_mat_t D(size,size);
//    gpu_mat_t E(size,size);

//     D = A+B;
//     viennacl::distributed::scheduler::finish();
//     E = A-B;
//     viennacl::distributed::scheduler::finish();
//     C = viennacl::generator::prod(D,E);
//     viennacl::distributed::scheduler::finish();

//    process_mem_usage(vm_usage,rss);
//    t.start();

//    D = A+B;
//    viennacl::distributed::scheduler::finish();
//    E = A-B;
//    viennacl::distributed::scheduler::finish();
//    C = viennacl::generator::prod(D,E);
//    viennacl::distributed::scheduler::finish();


//    return t.get();

     C=viennacl::generator::prod(A+B,A-B);
     viennacl::distributed::scheduler::finish();

     process_mem_usage(vm_usage,rss);
     t.start();

     C=viennacl::generator::prod(A+B,A-B);
     viennacl::distributed::scheduler::finish();

     return t.get();

}

int main(){
    viennacl::distributed::scheduler::add_all_available_devices(CL_DEVICE_TYPE_GPU);
    std::cout << "#Size \t Execution Time \t Memory Usage (MB) " << std::endl;
    for(unsigned int size = SIZE_INC ; size <= SIZE_INC * N_SIZES ; size+=SIZE_INC){
        double vm_usage;
        double rss;
        double time = run_benchmark<double>(size,vm_usage,rss);
        std::cout << size << "\t" << time << "\t" << rss/1024 << std::endl;
    }
}
