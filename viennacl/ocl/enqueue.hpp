#ifndef VIENNACL_OCL_ENQUEUE_HPP_
#define VIENNACL_OCL_ENQUEUE_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file enqueue.hpp
    @brief Enqueues kernels into command queues
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/command_queue.hpp"

namespace viennacl
{
  namespace generator{
      class custom_operation;
      void enqueue_custom_op(viennacl::generator::custom_operation & op, viennacl::ocl::command_queue const & queue);
  }
  
  namespace ocl
  { 
    template <typename KernelType>
    inline const event &enqueue(KernelType & k){
      return viennacl::ocl::find_context(viennacl::ocl::kernel::info<CL_KERNEL_CONTEXT>(k)).get_queue().enqueue(k);
    }

//    /** @brief Convenience function that enqueues the provided kernel into the first queue of the currently active device in the currently active context */
//    template <typename KernelType>
//    void enqueue(KernelType & k)
//    {
//      enqueue(k, viennacl::ocl::current_context().get_queue());
//    }
    
    inline void enqueue(viennacl::generator::custom_operation & op, viennacl::ocl::command_queue const & queue)
    {
      generator::enqueue_custom_op(op,queue);
    }

    inline void enqueue(viennacl::generator::custom_operation & op)
    {
      enqueue(op, viennacl::ocl::current_context().get_queue());
    }
    
  } // namespace ocl
} // namespace viennacl
#endif
