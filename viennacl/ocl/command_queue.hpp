#ifndef VIENNACL_OCL_COMMAND_QUEUE_HPP_
#define VIENNACL_OCL_COMMAND_QUEUE_HPP_

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

/** @file command_queue.hpp
    @brief Implementations of command queue representations
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <vector>
#include <string>
#include <sstream>
#include "viennacl/ocl/context.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/event.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief A class representing a command queue
    *
    */
    class command_queue
    {
      public:
        command_queue(viennacl::ocl::context * context) : context_(context){}
        command_queue(viennacl::ocl::context * context,viennacl::ocl::handle<cl_command_queue> h, cl_device_id dev) : context_(context), handle_(h) {}
        
//        //Copy constructor:
//        command_queue(command_queue const & other)
//        {
//          handle_ = other.handle_;
//          last_event_ = other.last_event_;
//        }

//        //assignment operator:
//        command_queue & operator=(command_queue const & other)
//        {
//          handle_ = other.handle_;
//          last_event_ = other.last_event_;
//          return *this;
//        }
        
        /** @brief Waits until all kernels in the queue have finished their execution */
        void finish() const
        {
          clFinish(handle_.get());
        }
        
        /** @brief Waits until all kernels in the queue have started their execution */
        void flush() const
        {
          clFlush(handle_.get());
        }

        viennacl::ocl::handle<cl_command_queue> const & handle() const { return handle_; }
        viennacl::ocl::handle<cl_command_queue>       & handle()       { return handle_; }

        /** @brief Creates an event for a specific kernel */
        void create_event(viennacl::ocl::kernel const & k)
        {
            cl_device_id id;
            cl_int err = clGetCommandQueueInfo(handle_.get(),  CL_QUEUE_DEVICE, sizeof(id), &id,NULL);
            VIENNACL_ERR_CHECK(err);
            if(events_.erase(k.handle().get())){
//                std::cout << "Removing previous event" << std::endl;
            }
            events_.insert(std::make_pair(k.handle().get(),viennacl::ocl::event(*this)));
        }

        /** @brief Deletes the event associated with a specific kernel  - assuming it exists */
        viennacl::ocl::event & get_event(viennacl::ocl::kernel const & k)
        {
            std::map<cl_kernel, viennacl::ocl::event>::iterator found = events_.find(k.handle().get());
            assert(found != events_.end() && "! No event for this kernel");
            return found->second;
        }

        /** @brief Returns the last event associated with an enqueueing */
        viennacl::ocl::event * last_event() { return last_event_; }

        /** @brief Enqueues a kernel in the provided queue */
        template <typename KernelType>
        viennacl::ocl::event & enqueue(KernelType & k)
        {

            create_event(k);
            viennacl::ocl::event & event = get_event(k);
            cl_event * event_ptr = const_cast<cl_event *>(&event.handle().get());
            last_event_ = &event;

            // 1D kernel:
            if (k.local_work_size(1) == 0)
            {
                #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
                std::cout << "ViennaCL: Starting 1D-kernel '" << k.name() << "'..." << std::endl;
                std::cout << "ViennaCL: Global work size: '"  << k.global_work_size() << "'..." << std::endl;
                std::cout << "ViennaCL: Local work size: '"   << k.local_work_size() << "'..." << std::endl;
                #endif

                size_t tmp_global = k.global_work_size();
                size_t tmp_local = k.local_work_size();

                cl_int err;
                if (tmp_global == 1 && tmp_local == 1)
                    err = clEnqueueTask(handle_.get(), k.handle().get(), 0, NULL, event_ptr);
                else
                    err = clEnqueueNDRangeKernel(handle_.get(), k.handle().get(), 1, NULL, &tmp_global, &tmp_local, 0, NULL, event_ptr);

                if (err != CL_SUCCESS)  //if not successful, try to start with smaller work size
                {
                    //std::cout << "FAIL: " << std::endl; exit(0);
                    while (err != CL_SUCCESS && tmp_local > 1)
                    {
                        //std::cout << "Flushing queue, then enqueuing again with half the size..." << std::endl;
                        //std::cout << "Error code: " << err << std::endl;

                        tmp_global /= 2;
                        tmp_local /= 2;

                        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
                        std::cout << "ViennaCL: Kernel start failed for '" << k.name() << "'." << std::endl;
                        std::cout << "ViennaCL: Global work size: '"  << tmp_global << "'..." << std::endl;
                        std::cout << "ViennaCL: Local work size: '"   << tmp_local << "'..." << std::endl;
                        #endif

                        finish();
                        err = clEnqueueNDRangeKernel(handle_.get(), k.handle().get(), 1, NULL, &tmp_global, &tmp_local, 0, NULL, event_ptr);
                    }

                    if (err != CL_SUCCESS)
                    {
                        //could not start kernel with any parameters
                        std::cerr << "ViennaCL: FATAL ERROR: Kernel start failed for '" << k.name() << "'." << std::endl;
                        std::cerr << "ViennaCL: Smaller work sizes could not solve the problem. " << std::endl;
                        VIENNACL_ERR_CHECK(err);
                    }
                    else
                    {
                        //remember parameters:
                        k.local_work_size(0, tmp_local);
                        k.global_work_size(0, tmp_global);
                        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
                        std::cout << "ViennaCL: Kernel '" << k.name() << "' now uses global work size " << tmp_global << " and local work size " << tmp_local << "."  << std::endl;
                        #endif
                    }
                }
            }
            else //2D kernel
            {
                #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
                std::cout << "ViennaCL: Starting 2D-kernel '" << k.name() << "'..." << std::endl;
                std::cout << "ViennaCL: Global work size: '"  << k.global_work_size(0) << ", " << k.global_work_size(1) << "'..." << std::endl;
                std::cout << "ViennaCL: Local work size: '"   << k.local_work_size(0) << ", " << k.local_work_size(1) << "'..." << std::endl;
                #endif

                size_t tmp_global[2];
                tmp_global[0] = k.global_work_size(0);
                tmp_global[1] = k.global_work_size(1);

                size_t tmp_local[2];
                tmp_local[0] = k.local_work_size(0);
                tmp_local[1] = k.local_work_size(1);

                cl_int err = clEnqueueNDRangeKernel(handle_.get(), k.handle().get(), 2, NULL, tmp_global, tmp_local, 0, NULL, event_ptr);

                if (err != CL_SUCCESS)
                {
                    //could not start kernel with any parameters
                    std::cerr << "ViennaCL: FATAL ERROR: Kernel start failed for '" << k.name() << "'." << std::endl;
                    VIENNACL_ERR_CHECK(err);
                }

          }

          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
          finish();
          std::cout << "ViennaCL: Kernel " << k.name() << " finished!" << std::endl;
          #endif

          return event;

        } //enqueue()
      private:
        
        viennacl::ocl::handle<cl_command_queue> handle_;

        viennacl::ocl::context * context_;
        viennacl::ocl::event * last_event_;
        std::map<cl_kernel, viennacl::ocl::event> events_;
    };

 
    
  } //namespace ocl
} //namespace viennacl

#endif
