#ifndef VIENNACL_DISTRIBUTED_UTILS_HPP_
#define VIENNACL_DISTRIBUTED_UTILS_HPP_

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

/** @file utils.hpp
    @brief Implementation of several utils
*/

#include "viennacl/tools/tools.hpp"
#include "boost/shared_ptr.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/ocl/backend.hpp"

namespace viennacl{

namespace distributed{

namespace utils{

template <class INT_TYPE>
INT_TYPE roundDownToPreviousMultiple(INT_TYPE to_reach, INT_TYPE base)
{
  if (to_reach % base == 0) return to_reach;
  return (to_reach / base) * base;
}

template<class ScalarType>
vcl_size_t matrix_block_size(){
    vcl_size_t mem_chunk_size = 128*1024*1024;
    return roundDownToPreviousMultiple<vcl_size_t>(sqrt(mem_chunk_size/(sizeof(ScalarType))),64);
}

template<class ScalarType>
vcl_size_t vector_block_size(){
    return matrix_block_size<ScalarType>()*matrix_block_size<ScalarType>();
}


/** @brief Conversion between viennacl layout types and boost layout types */
template<class T>
struct layout_wrapper;

template<>
struct layout_wrapper<viennacl::row_major>{ typedef boost::numeric::ublas::row_major Result ; };

template<>
struct layout_wrapper<viennacl::column_major>{ typedef boost::numeric::ublas::column_major Result ; };

template<class T>
struct get_cpu_type;

template<class SCALARTYPE, class F, unsigned int Alignment>
struct get_cpu_type<viennacl::matrix<SCALARTYPE,F,Alignment> >{
    typedef boost::numeric::ublas::matrix<SCALARTYPE, typename layout_wrapper<F>::Result> type;
};

template<class SCALARTYPE, class F, unsigned int Alignment>
struct get_cpu_type<const viennacl::matrix<SCALARTYPE,F,Alignment> >{
    typedef const boost::numeric::ublas::matrix<SCALARTYPE, typename layout_wrapper<F>::Result> type;
};


template<class SCALARTYPE, unsigned int Alignment>
struct get_cpu_type<viennacl::vector<SCALARTYPE,Alignment> >{
    typedef std::vector<SCALARTYPE> type;
};

template<class SCALARTYPE,unsigned int Alignment>
struct get_cpu_type<const viennacl::vector<SCALARTYPE,Alignment> >{
    typedef const std::vector<SCALARTYPE> type;
};

/** @brief Storage for lazy gpu allocation */

template<class T, class CpuT>
struct alloc_impl;

template<class ScalarType, class F, class CpuT>
struct alloc_impl<viennacl::matrix<ScalarType, F>, CpuT>{
    viennacl::matrix<ScalarType, F>* operator()(CpuT const & cpu_data){
        size_t size1=cpu_data.size1(), size2=cpu_data.size2();
        viennacl::matrix<ScalarType, F>* p = new viennacl::matrix<ScalarType, F>(size1, size2);
        cl_mem h = p->handle().opencl_handle();
        clEnqueueWriteBuffer(viennacl::ocl::current_context().get_queue().handle().get(),h,true,0,size1*size2,&cpu_data(0,0),0,NULL,NULL);
        return p;
    }
};

template<class ScalarType, class F, class CpuT>
struct alloc_impl<const viennacl::matrix<ScalarType, F>, CpuT>{
    const viennacl::matrix<ScalarType, F>* operator()(CpuT const & cpu_data){
        size_t size1=cpu_data.size1(), size2=cpu_data.size2();
        const viennacl::matrix<ScalarType, F>* p = new viennacl::matrix<ScalarType, F>(size1, size2);
        cl_mem h = p->handle().opencl_handle();
        clEnqueueWriteBuffer(viennacl::ocl::current_context().get_queue().handle().get(),h,true,0,size1*size2,&cpu_data(0,0),0,NULL,NULL);
        return p;
    }
};

//template<class ScalarType, class F, class CpuT>
//struct alloc_impl<const viennacl::matrix<ScalarType, F>, CpuT>{
//    viennacl::matrix<ScalarType, F>* operator()(viennacl::ocl::context &ctxt, CpuT const & cpu_data){
//        return new viennacl::matrix<ScalarType, F>(cpu_data.size1(), cpu_data.size2(), &cpu_data(0,0), ctxt);
//    }
//};


//template<class ScalarType, class F, class CpuT>
//viennacl::matrix<ScalarType, F>* alloc_impl<viennacl::matrix<ScalarType, F>, CpuT>(viennacl::ocl::context & ctxt, CpuT const & cpu_data){
//    return new viennacl::matrix<ScalarType, F>(cpu_data.size1(), cpu_data.size2(), &cpu_data(0,0), ctxt);
//}

//template<class ScalarType>
//viennacl::vector<ScalarType>* alloc_impl(viennacl::ocl::context & ctxt){
//    return new viennacl::vector<ScalarType>(cpu_data.size1(), cpu_data.size2(), &cpu_data(0,0), ctxt);
//}

template<class T>
class gpu_wrapper{
public:
    typedef typename get_cpu_type<T>::type cpu_t;
    gpu_wrapper(cpu_t & _cpu_data) : cpu_data(_cpu_data){ }

    gpu_wrapper(gpu_wrapper const & other) : cpu_data(other.cpu_data){
        assert(gpu_structure_.get() == NULL);
    }

    void alloc(){
        gpu_structure_.reset(alloc_impl<T,cpu_t>()(cpu_data));
    }

    void free(){
        gpu_structure_.reset();
    }

    void transfer_back(){
        size_t internal_size = cpu_data.size1() * cpu_data.size2();
        clEnqueueReadBuffer(viennacl::ocl::current_context().get_queue().handle().get(),gpu_structure_->handle().opencl_handle(),true,0,internal_size,&cpu_data(0,0),0,NULL,NULL);
    }

    T * gpu_structure_ptr(){
        return gpu_structure_.get();
    }


private:
    cpu_t & cpu_data;
    boost::shared_ptr<T> gpu_structure_;
};


}

}

}

#endif // VIENNACL_DISTRIBUTED_UTILS_HPP_
