#ifndef VIENNACL_DISTRIBUTED_TASK_HPP_
#define VIENNACL_DISTRIBUTED_TASK_HPP_


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

#include "viennacl/distributed/utils.hpp"
#include "boost/thread/mutex.hpp"
#include <boost/type_traits/is_const.hpp>
#include "viennacl/ocl/utils.hpp"
#include "viennacl/ocl/event.hpp"

/** @file task.hpp
    @brief Implementation of a task
*/

namespace viennacl{

namespace distributed{

template<class GPU_WRAPPER_T>
class transfer_handler;

template<class T>
class transfer_handler< viennacl::distributed::utils::gpu_wrapper<T> >{
private:
    template<class MAT_T>
    void transfer_back(viennacl::distributed::utils::gpu_wrapper<MAT_T const> & ){ }

    template<class MAT_T>
    void transfer_back(viennacl::distributed::utils::gpu_wrapper<MAT_T>  & ){
        wrapper_.transfer_back();
    }
public:
    transfer_handler(viennacl::distributed::utils::gpu_wrapper<T> & wrapper) : wrapper_(wrapper){
        wrapper_.alloc();
    }

    T & gpu_structure(){
        return *wrapper_.gpu_structure_ptr();
    }

    ~transfer_handler(){
        transfer_back(wrapper_);
        wrapper_.free();
    }

private:
    viennacl::distributed::utils::gpu_wrapper<T> & wrapper_;
};


class task{
protected:
public:
    virtual viennacl::ocl::event * run() = 0;

    std::string const & info() const{
        return info_;
    }

    void info(std::string const & new_info){
        info_ = new_info;
    }

    virtual ~task(){    }

protected:
    std::string info_;
};


template<class ARG0, class ARG1, class RES>
class task2 : public task{
private:
    typedef std::function<void (RES&, ARG0 const &, ARG1 const &)> fun_t;
    typedef viennacl::distributed::utils::gpu_wrapper<const ARG0> gwrap0_t;
    typedef viennacl::distributed::utils::gpu_wrapper<const ARG1> gwrap1_t;
    typedef viennacl::distributed::utils::gpu_wrapper<RES> gwrapres_t;
public:
    task2(fun_t fun, gwrap0_t arg0, gwrap1_t arg1, gwrapres_t  res) : arg0_(arg0),arg1_(arg1), res_(res), fun_(fun){ }

    viennacl::ocl::event * run(){
#ifdef VIENNACL_DEBUG_SCHEDULER
        std::cout << "Running " << info() << std::endl;
#endif
        transfer_handler< gwrap0_t > gpu_arg0(arg0_);
        transfer_handler< gwrap1_t > gpu_arg1(arg1_);
        transfer_handler< gwrapres_t > gpu_res(res_);

        Timer t;
        t.start();

        fun_(gpu_res.gpu_structure(), gpu_arg0.gpu_structure(), gpu_arg1.gpu_structure());
        return viennacl::ocl::get_queue().last_event();
    }
private:
    gwrap0_t arg0_;
    gwrap1_t arg1_;
    gwrapres_t  res_;
    fun_t fun_;
};

}

}
#endif // TASK_HPP
