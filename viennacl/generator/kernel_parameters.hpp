#ifndef VIENNACL_GENERATOR_KERNEL_PARAMETERS_HPP
#define VIENNACL_GENERATOR_KERNEL_PARAMETERS_HPP

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


/** @file kernel_parameters.hpp
    @brief This file holds the code necessary for reading kernel parameters from XML files using pugixml
*/

#include "viennacl/ocl/backend.hpp"
#include "pugixml/src/pugixml.hpp"

namespace viennacl
{
  namespace generator
  {
    namespace tag
    {
      static std::string root     = "parameters";
      static std::string devices  = "devices";
      static std::string device   = "device";
      static std::string name     = "name";
      static std::string driver   = "driver";
      static std::string numeric  = "numeric";
      static std::string alignment = "alignment";
    } // end namespace tag

    namespace val {
      static std::string globsize = "globalsize";
      static std::string locsize  = "localsize";
      static std::string vec      = "vector";
      static std::string matrix   = "matrix";
      static std::string fl       = "float";
      static std::string dbl      = "double";
    }

    namespace blas3_tag{
      static std::string ml       = "ml";
      static std::string kl       = "kl";
      static std::string nl       = "nl";
      static std::string ms       = "ms";
      static std::string ks       = "ks";
      static std::string ns       = "ns";
      static std::string layout   = "layout";
      static std::string transposed = "transposed";
    }

    namespace blas3_val{
      static std::string row_major   = "row_major";
      static std::string column_major   = "column_major";
    }



  } // end namespace generator

} // end namespace viennacl

#endif

