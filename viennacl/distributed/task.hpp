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
#include "viennacl/generator/forwards.h"

/** @file task.hpp
    @brief Implementation of a task
*/

namespace viennacl{

namespace distributed{

//template<class GPU_WRAPPER_T>
//class transfer_handler;

//template<class T>
//class transfer_handler< viennacl::distributed::utils::gpu_wrapper<T> >{
//private:
//    template<class MAT_T>
//    void transfer_back(viennacl::distributed::utils::gpu_wrapper<MAT_T const> & ){ }

//    template<class MAT_T>
//    void transfer_back(viennacl::distributed::utils::gpu_wrapper<MAT_T>  & ){
//        wrapper_.transfer_back();
//    }
//public:
//    transfer_handler(viennacl::distributed::utils::gpu_wrapper<T> & wrapper) : wrapper_(wrapper){
//        wrapper_.alloc();
//    }

//    T & gpu_structure(){
//        return *wrapper_.gpu_structure_ptr();
//    }

//    ~transfer_handler(){
//        transfer_back(wrapper_);
//        wrapper_.free();
//    }

//private:
//    viennacl::distributed::utils::gpu_wrapper<T> & wrapper_;
//};

template<class T>
struct wrapper_to_matrix;

template<class T>
struct wrapper_to_matrix<utils::gpu_wrapper<T> >{
    typedef generator::dummy_matrix<T> type;
};

class wrapper_to_matrix_fun{
public:
    template<class T>
    generator::dummy_matrix<T> operator()(utils::gpu_wrapper<T> const & t) const{
        return generator::dummy_matrix<T>(*t.gpu_structure_ptr());
    }
};

struct alloc_fun{
    template<class T>
    void operator()(utils::gpu_wrapper<T> const & t){ t.alloc(); }
};

struct free_fun{
    template<class T>
    void operator()(utils::gpu_wrapper<T> const & t){ t.free(); }
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




template<class T, class A>
class task2 : public task{
private:
    typedef std::function<T> fun_t;
    typedef typename utils::make_deep_copy<A>::result_type ArgsT;
public:
    task2(fun_t fun, A const & args) : fun_(fun), args_(args){ }

    viennacl::ocl::event * run(){
#ifdef VIENNACL_DEBUG_SCHEDULER
        std::cout << "Running " << info() << std::endl;
#endif
        utils::execute<A>()(args_,alloc_fun());
        typedef typename utils::transform<A,wrapper_to_matrix,wrapper_to_matrix_fun>::result_type Tmp2Type;
        Tmp2Type tmp2(utils::transform<A,wrapper_to_matrix,wrapper_to_matrix_fun>::result(args_,wrapper_to_matrix_fun()));

        fun_(utils::make_shallow_copy<Tmp2Type>::result(tmp2));

        args_.lhs().transfer_back();
        utils::execute<A>()(args_,free_fun());
        return viennacl::ocl::get_queue().last_event();
    }
private:
    fun_t fun_;
    A args_;
};

}

}
#endif // TASK_HPP
