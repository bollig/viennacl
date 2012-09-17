#ifndef VIENNACL_LINALG_PARALLEL_NORM_INF_HPP_
#define VIENNACL_LINALG_PARALLEL_NORM_INF_HPP_

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

/** @file norm_inf.hpp
    @brief Generic interface for the l^inf-norm. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include <math.h>    //for sqrt()
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"

#include "viennacl/linalg/norm_inf.hpp"
#include "utils/comm/communicator.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/filesystem.hpp>

namespace viennacl
{
  //
  // generic norm_inf function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {

    // Generic for GPU (compute norm on gpu, transfer scalar to cpu, reduce and return to the gpu)
    template <typename VectorType>
    typename VectorType::value_type
    parallel_norm_inf(VectorType const& v, Communicator const& comm_unit, typename VectorType::value_type& dummy) 
    {
        // Force the GPU norm to the cpu. 
        double norm = viennacl::linalg::norm_inf(v);
        double r[1]; 
        
        MPI_Allreduce(&norm, r, 1, MPI_DOUBLE, MPI_MAX, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 
       
        return (typename VectorType::value_type)r[0]; 
    }


    template <typename VectorType>
    //typename VectorType::value_type
    double
    parallel_norm_inf(VectorType const& v, Communicator const& comm_unit, double& dummy) 
    {
        // MPI will segfault if we dont declare norm as an array type
        //typename VectorType::value_type norm = viennacl::linalg::norm_inf(v);
        double norm = viennacl::linalg::norm_inf(v);
        double r[1]; 
        
        MPI_Allreduce(&norm, r, 1, MPI_DOUBLE, MPI_MAX, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 
       
        return r[0]; 
    }

    template <typename VectorType>
    typename VectorType::value_type
    norm_inf(VectorType const& v, Communicator const& comm_unit) 
    {
        // Redirect the call to appropriate float or double routines (Type size matters with MPI)
        typename VectorType::value_type dummy; 
        return parallel_norm_inf(v, comm_unit, dummy ); 
    }

    template < >
    double
    norm_inf(boost::numeric::ublas::vector_range< boost::numeric::ublas::vector<double> > const& v, Communicator const& comm_unit) 
    {
            double norm = boost::numeric::ublas::norm_inf(v);
            double r[1]; 

            MPI_Allreduce(&norm, r, 1, MPI_DOUBLE, MPI_MAX, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

            return r[0]; 
    }

    template < >
    double
    norm_inf(boost::numeric::ublas::vector_slice< boost::numeric::ublas::vector<double> > const& v, Communicator const& comm_unit) 
    {
            double norm = boost::numeric::ublas::norm_inf(v);
            double r[1]; 

            MPI_Allreduce(&norm, r, 1, MPI_DOUBLE, MPI_MAX, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

            return r[0]; 
    }
 } // end namespace linalg
} // end namespace viennacl
#endif





