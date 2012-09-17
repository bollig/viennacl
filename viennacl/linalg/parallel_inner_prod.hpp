#ifndef VIENNACL_LINALG_PARALLEL_INNER_PROD_HPP_
#define VIENNACL_LINALG_PARALLEL_INNER_PROD_HPP_

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

/** @file inner_prod.hpp
  @brief Generic interface for the computation of inner products. See viennacl/linalg/vector_operations.hpp for implementations.
  */

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"

#include "viennacl/linalg/inner_prod.hpp"
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
    // generic inner_prod function
    //   uses tag dispatch to identify which algorithm
    //   should be called 
    //
    namespace linalg 
    {
        // Generic for GPU (compute norm on gpu, transfer scalar to cpu, reduce and return to the gpu)
        template <typename VectorType, typename VectorType2>
            typename VectorType::value_type
            parallel_inner_prod(VectorType const& v, VectorType2 const& v2, Communicator const& comm_unit, typename VectorType::value_type& dummy) 
            {
                // Force the GPU norm to the cpu. 
                double val = viennacl::linalg::inner_prod(v, v2);
                double r[1]; 

                MPI_Allreduce(&val, r, 1, MPI_DOUBLE, MPI_SUM, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

                return (typename VectorType::value_type)r[0]; 

            }


#if 0
        template <typename VectorType, typename VectorType2>
            //typename VectorType::value_type
            double
            parallel_inner_prod(VectorType const& v, VectorType2 const& v2, Communicator const& comm_unit, double& dummy) 
            {
                double val = viennacl::linalg::inner_prod(v, v2);
                double r[1]; 

                MPI_Allreduce(&val, r, 1, MPI_DOUBLE, MPI_SUM, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

                return r[0]; 
            }
#endif 
 

        template < typename VectorType , typename VectorType2 >
            double 
            parallel_inner_prod( VectorType const& v, viennacl::vector_range< VectorType2 > const& v2, Communicator const& comm_unit, typename VectorType::value_type & dummy) 
            {
                // Force the GPU norm to the cpu. 
                double val = viennacl::linalg::inner_prod(v, v2);
                double r[1]; 

                MPI_Allreduce(&val, r, 1, MPI_DOUBLE, MPI_SUM, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

                return r[0]; 

            }


        template < typename VectorType, typename VectorType2 >
            double 
            parallel_inner_prod( VectorType const& v, boost::numeric::ublas::vector_range< VectorType2 > const& v2, Communicator const& comm_unit, typename VectorType::value_type & dummy) 
            {
                // Force the GPU norm to the cpu. 
                double val = boost::numeric::ublas::inner_prod(v, v2);
                double r[1]; 

                MPI_Allreduce(&val, r, 1, MPI_DOUBLE, MPI_SUM, comm_unit.getComm()); 

                MPI_Barrier(MPI_COMM_WORLD); 
                return r[0]; 

            }



#if 0
        template <typename VectorType, typename VectorType2>
            float 
            parallel_inner_prod(VectorType const& v, VectorType2 const& v2, Communicator const& comm_unit, float& dummy) 
            {
                float val = viennacl::linalg::inner_prod(v, v2);
                float r[1]; 

                MPI_Allreduce(&val, r, 1, MPI_FLOAT, MPI_SUM, comm_unit.getComm()); 
            MPI_Barrier(MPI_COMM_WORLD); 

                return r[0]; 
            }
#endif 

        template <typename VectorType, typename VectorType2>
            typename VectorType::value_type
            inner_prod(VectorType const& v, VectorType2 const& v2, Communicator const& comm_unit) 
            {
                // Redirect the call to appropriate float or double routines (Type size matters with MPI)
                typename VectorType::value_type dummy; 
                return parallel_inner_prod(v, v2, comm_unit, dummy ); 
            }

#if 0
        template< typename VectorT1, typename VectorT2 >
            typename VectorT1::value_type
            inner_prod(VectorT1 const& v1, VectorT2 const& v2, 
                    typename viennacl::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                    >::type* dummy = 0)
            {
                //std::cout << "ublas .. " << std::endl;
                return boost::numeric::ublas::inner_prod(v1, v2);
            }
#endif 
    } // end namespace linalg
} // end namespace viennacl
#endif


