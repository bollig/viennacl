#ifndef VIENNACL_OCL_BACKEND_HPP_
#define VIENNACL_OCL_BACKEND_HPP_

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

/** @file backend.hpp
    @brief Implementations of the OpenCL backend, where all contexts are stored in.
*/

#include <vector>
#include "viennacl/ocl/context.hpp"
#include "viennacl/ocl/enqueue.hpp"

#if __cplusplus > 199711L
#include <thread>
#endif

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief A backend that provides contexts for ViennaCL objects (vector, matrix, etc.) */
    template <bool dummy = false>  //never use parameter other than default (introduced for linkage issues only)
    class backend
    {
      public:
        /** @brief Switches the current context to the context identified by i 
        *
        * @param i   ID of the new active context
        */

#if __cplusplus > 199711L
        static void switch_context(long i){
            current_context_id_[std::this_thread::get_id()] = i;
        }
#else
        static void switch_context(long i)
        {
          current_context_id_ = i;
        }
#endif

        
        /** @brief Returns the current active context */
#if __cplusplus > 199711L
        static viennacl::ocl::context & current_context(std::thread::id id = std::this_thread::get_id()){
          long context_id = current_context_id_[id];
#else
        static viennacl::ocl::context & current_context(){
          long context_id = current_context_id_;
#endif

          if (!initialized_[context_id])
          {
            //std::cout << "Initializing context no. " << context_id << std::endl;
            contexts_[context_id].init();
            //create one queue per device:
            std::vector<viennacl::ocl::device> devices = contexts_[context_id].devices();
            for (size_t j = 0; j<devices.size(); ++j)
              contexts_[context_id].add_queue(devices[j]);
            initialized_[context_id] = true;
            /*
            std::cout << "Context no. " << context_id << " initialized with " << devices.size() << " devices" << std::endl;
            std::cout << "Device id: " << devices[0].id() << std::endl;
            std::cout << "Current device id: " << contexts_[context_id].current_device().id() << std::endl; */
          }
          return contexts_[context_id];
        }
        
        /** @brief Returns the current queue for the active device in the active context */
        static viennacl::ocl::command_queue & get_queue()
        {
          return current_context().get_queue();
        }

        /** @brief Sets a number of devices for the context.
        *
        * @param i    ID of the context to be set up
        * @param devices A vector of OpenCL device-IDs that should be added to the context
        */
        static void setup_context(long i,
                                  std::vector<cl_device_id> const & devices)
        {
          if (initialized_[i])
            std::cerr << "ViennaCL: Warning in init_context(): Providing a list of devices has no effect, because context for ViennaCL is already created!" << std::endl;
          else
          {
            //set devices for context:
            for (size_t j = 0; j<devices.size(); ++j)
              contexts_[i].add_device(devices[j]);
          }
        }

        /** @brief Initializes ViennaCL with an already existing context
        *
        * @param i    ID of the context to be set up
        * @param c    The OpenCL handle of the existing context
        * @param devices A vector of OpenCL device-IDs that should be added to the context
        * @param queues   A map of queues for each device
        */
        static void setup_context(long i,
                                  cl_context c,
                                  std::vector<cl_device_id> const & devices,
                                  std::map< cl_device_id, std::vector< cl_command_queue > > const & queues)
        {
          assert(devices.size() == queues.size() && bool("ViennaCL expects one queue per device!"));
          
          if (initialized_[i])
            std::cerr << "ViennaCL: Warning in init_context(): Providing a list of devices has no effect, because context for ViennaCL is already created!" << std::endl;
          else
          {
            //set devices for context:
            for (size_t j = 0; j<devices.size(); ++j)
              contexts_[i].add_device(devices[j]);
            
            //init context:
            contexts_[i].init(c);
            
            //add queues:
            typedef typename std::map< cl_device_id, std::vector< cl_command_queue > >::const_iterator queue_iterator;
            for (queue_iterator qit = queues.begin();
                                qit != queues.end();
                              ++qit)
            {
              std::vector<cl_command_queue> const & queues_for_device = qit->second;
              for (size_t j=0; j<queues_for_device.size(); ++j)
                contexts_[i].add_queue(qit->first, queues_for_device[j]);
            }
            
            initialized_[i] = true;
          }
        }

        /** @brief Initializes ViennaCL with an already existing context
        *
        * @param i    ID of the context to be set up
        * @param c    The OpenCL handle of the existing context
        * @param devices A vector of OpenCL device-IDs that should be added to the context
        * @param queue   One queue per device
        */
        static void setup_context(long i, cl_context c, std::vector<cl_device_id> const & devices, std::vector<cl_command_queue> const & queue)
        {
          assert(devices.size() == queue.size() && bool("ViennaCL expects one queue per device!"));
          
          //wrap queue vector into map
          std::map< cl_device_id, std::vector<cl_command_queue> > queues_map;
          for (size_t j = 0; j<devices.size(); ++j)
            queues_map[devices[j]].push_back(queue[j]);
          
          setup_context(i, c, devices, queues_map);
        }

        /** @brief Sets the context device type */
        static void set_context_device_type(long i, cl_device_type t)
        {
          contexts_[i].default_device_type(t);
        }

        /** @brief Sets the context device type */
        static void set_context_platform_index(long i, std::size_t pf_index)
        {
          contexts_[i].platform_index(pf_index);
        }

        static viennacl::ocl::context& find_context(cl_context handle){
            for(std::map<long,viennacl::ocl::context>::iterator it = contexts_.begin(); it!=contexts_.end();++it){
                if(it->second.handle().get()==handle) return it->second;
            }
            throw "No such context in the backend";
        }



      private:
#if __cplusplus > 199711L
        static std::map<std::thread::id, long> current_context_id_;
#else
        static long current_context_id_;
#endif
        static std::map<long, bool> initialized_;
        static std::map<long, viennacl::ocl::context> contexts_;
    };
    

#if __cplusplus > 199711L
//    static std::map<std::thread::id, long> initialize_context_id(){
//        std::map<std::thread::id, long> res;
//        res.insert(std::make_pair(std::this_thread::get_id(),0));
//        return res;
//    }

    template <bool dummy>
    std::map<std::thread::id, long> backend<dummy>::current_context_id_;
#else
    template <bool dummy>
    long backend<dummy>::current_context_id_ = 0;
#endif

    template <bool dummy>
    std::map<long, bool> backend<dummy>::initialized_;

    template <bool dummy>
    std::map<long, viennacl::ocl::context> backend<dummy>::contexts_;
    
    ////////////////////// current context //////////////////
    /** @brief Convenience function for returning the current context */
    inline viennacl::ocl::context & current_context()
    {
      return viennacl::ocl::backend<>::current_context();
    }

    /** @brief Convenience function for switching the current context */
    inline void switch_context(long i)
    {
      viennacl::ocl::backend<>::switch_context(i);
    }
    

    /** @brief Convenience function for setting devices for a context */
    inline void setup_context(long i,
                              std::vector<cl_device_id> const & devices)
    {
      viennacl::ocl::backend<>::setup_context(i, devices);
    }

    /** @brief Convenience function for setting up a context in ViennaCL from an existing OpenCL context */
    inline void setup_context(long i,
                              cl_context c,
                              std::vector<cl_device_id> const & devices,
                              std::map< cl_device_id, std::vector<cl_command_queue> > const & queues)
    {
      viennacl::ocl::backend<>::setup_context(i, c, devices, queues);
    }
    
    /** @brief Convenience function for setting up a context in ViennaCL from an existing OpenCL context */
    inline void setup_context(long i, cl_context c, std::vector<cl_device_id> const & devices, std::vector<cl_command_queue> const & queues)
    {
      viennacl::ocl::backend<>::setup_context(i, c, devices, queues);
    }

    /** @brief Convenience function for setting up a context in ViennaCL from an existing OpenCL context */
    inline void setup_context(long i, cl_context c, cl_device_id d, cl_command_queue q)
    {
      std::vector<cl_device_id> devices(1);
      std::vector<cl_command_queue> queues(1);
      devices[0] = d;
      queues[0] = q;
      viennacl::ocl::backend<>::setup_context(i, c, devices, queues);
    }

    /** @brief Convenience function for setting the default device type for a context */
    inline void set_context_device_type(long i, cl_device_type dev_type)
    {
      viennacl::ocl::backend<>::set_context_device_type(i, dev_type);
    }

    /** @brief Convenience function for setting the default device type for a context to GPUs */
    inline void set_context_device_type(long i, viennacl::ocl::gpu_tag)
    {
      set_context_device_type(i, CL_DEVICE_TYPE_GPU);
    }

    /** @brief Convenience function for setting the default device type for a context to CPUs */
    inline void set_context_device_type(long i, viennacl::ocl::cpu_tag)
    {
      set_context_device_type(i, CL_DEVICE_TYPE_CPU);
    }

    /** @brief Convenience function for setting the default device type for a context to the default OpenCL device type */
    inline void set_context_device_type(long i, viennacl::ocl::default_tag)
    {
      set_context_device_type(i, CL_DEVICE_TYPE_DEFAULT);
    }

    /** @brief Convenience function for setting the default device type for a context to accelerators */
    inline void set_context_device_type(long i, viennacl::ocl::accelerator_tag)
    {
      set_context_device_type(i, CL_DEVICE_TYPE_ACCELERATOR);
    }

    
    /** @brief Convenience function for setting the platform index
     * 
     * @param i         Context ID
     * @param pf_index  The platform index as returned by clGetPlatformIDs(). This is not the ID of type cl_platform_id!
     */
    inline void set_context_platform_index(long i, std::size_t pf_index)
    {
      viennacl::ocl::backend<>::set_context_platform_index(i, pf_index);
    }
    

    ///////////////////////// get queues ///////////////////
    /** @brief Convenience function for getting the default queue for the currently active device in the active context */
    inline viennacl::ocl::command_queue & get_queue()
    {
      return viennacl::ocl::current_context().get_queue();
    }
    
    /** @brief Convenience function for getting the queue for a particular device in the current active context */
    inline viennacl::ocl::command_queue & get_queue(viennacl::ocl::device d, unsigned int queue_id = 0)
    {
      return viennacl::ocl::current_context().get_queue(d.id(), queue_id);
    }

    /** @brief Convenience function for getting the queue for a particular device in the current active context */
    inline viennacl::ocl::command_queue & get_queue(cl_device_id dev_id, unsigned int queue_id = 0)
    {
      return viennacl::ocl::current_context().get_queue(dev_id, queue_id);
    }


    /** @brief Convenience function for getting the kernel for a particular program from the current active context */
    inline viennacl::ocl::kernel & get_kernel(std::string const & prog_name, std::string const & kernel_name)
    {
      return viennacl::ocl::current_context().get_program(prog_name).get_kernel(kernel_name);
    }

    /** @brief Convenience function for switching the active device in the current context */
    inline void switch_device(viennacl::ocl::device & d)
    {
      viennacl::ocl::current_context().switch_device(d);
    }
    
    /** @brief Convenience function for returning the active device in the current context */
    inline viennacl::ocl::device const & current_device()
    {
      return viennacl::ocl::current_context().current_device();
    }

    inline viennacl::ocl::context & find_context(cl_context handle){
        return viennacl::ocl::backend<>::find_context(handle);
    }

  } //ocl
} //viennacl
#endif
