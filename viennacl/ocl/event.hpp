#ifndef VIENNACL_OCL_EVENT_HPP_
#define VIENNACL_OCL_EVENT_HPP_

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

/** @file event.hpp
    @brief Implementations of event representations
*/

#include "CL/cl.h"

#include "list"

#include "viennacl/ocl/forwards.h"
#include "viennacl/ocl/handle.hpp"

namespace viennacl{

namespace ocl{

class event{
public:
    event(viennacl::ocl::command_queue & queue) : queue_(queue){ }

    cl_int status() const
    {
        cl_int res = 0;
        cl_int err = clGetEventInfo(handle_.get(),CL_EVENT_COMMAND_EXECUTION_STATUS,sizeof(res),&res,NULL);
        VIENNACL_ERR_CHECK(err);
        return res;
    }

    void init_callback(){
        cl_int err = clSetEventCallback(handle_.get(), CL_COMPLETE, callback_fn_, callback_data_);
        VIENNACL_ERR_CHECK(err);
    }

    void set_callback( void (CL_CALLBACK  *new_fn) (cl_event , cl_int, void *)){
        callback_fn_ = new_fn;
    }

    void set_callback_data( void* new_data){
        callback_data_ = new_data;
    }

    viennacl::ocl::handle<cl_event> const & handle() const{
        return handle_;
    }

    viennacl::ocl::handle<cl_event> & handle(){
        return handle_;
    }

private:
    viennacl::ocl::command_queue & queue_;
    viennacl::ocl::handle<cl_event> handle_;
    void (CL_CALLBACK  *callback_fn_) (cl_event , cl_int, void *);
    void *callback_data_;
};

#ifdef VIENNACL_PROFILING_ENABLE

cl_ulong time_of_execution_us(cl_event e){
    cl_ulong time_start, time_end;
    double total_time;
    clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
    clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
    total_time = time_end - time_start;
    return total_time/1000;
}

cl_ulong time_of_execution_us(viennacl::ocl::event & e){
    return time_of_execution_us(e.handle());
}

#endif

//static void wait_for_events(std::list<viennacl::ocl::event> const & data)
//{
//    cl_uint length = static_cast<cl_uint>(data.size());
//    cl_event * raw_data = new cl_event[length];
//    for(std::list<viennacl::ocl::event>::const_iterator it = data.begin() ; it != data.end() ; ++it){
//        unsigned int index = std::distance(it, data.begin());
//        raw_data[index] = it->handle().get();
//    }
//    cl_int err = clWaitForEvents(length,raw_data);
//    //Deletes before event checking to avoid memory leaks.
//    delete[] raw_data;
//    VIENNACL_ERR_CHECK(err);
//}

}

}
#endif // EVENT_HPP
