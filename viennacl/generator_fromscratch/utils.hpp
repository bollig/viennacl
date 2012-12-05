#ifndef VIENNACL_GENERATOR_UTILS_HPP
#define VIENNACL_GENERATOR_UTILS_HPP

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

#include <sstream>

namespace viennacl{
	
	namespace generator{
		
            template <class T>
            inline std::string to_string ( T const t )
            {
              std::stringstream ss;
              ss << t;
              return ss.str();
            }

            template<class T, class U>
            static bool is_type(U* p){
                return dynamic_cast<T *>(p);
            }


            template<class Base,class Target>
            struct Base2Target { Target* operator ()( Base* value ) const { return dynamic_cast<Target*>(value); } };


            template<class Base,class Target>
            struct UnsafeBase2Target { Target* operator ()( Base* value ) const { return static_cast<Target*>(value); } };


		    template<class T>
			struct print_align1_type;

			template<>
			struct print_align1_type<int> 
			{
			  static const std::string value() { return "int"; }
			};

			template<>
			struct print_align1_type<unsigned int>
			{
			  static const std::string value() { return "unsigned int"; }
			};

			template<>
			struct print_align1_type<long> 
			{
			  static const std::string value() { return "long"; }
			};

			template<>
			struct print_align1_type<unsigned long> 
			{
			  static const std::string value() { return "long"; }
			};

			template<>
			struct print_align1_type<float> 
			{
			  static const std::string value() { return "float"; }
			};

			template<>
			struct print_align1_type<double> 
			{
			  static const std::string value() { return "double"; }
			};

			template<typename T, unsigned int ALIGNMENT>
			struct print_aligned_type 
			{
				static const std::string value() 
				{
				return print_align1_type<T>::value() + to_string ( ALIGNMENT );
			  }
			};

			template<typename T>
			struct print_aligned_type<T, 1>
			{
				static const std::string value() 
				{
				return print_align1_type<T>::value();
			  }
			};

			template<typename T, unsigned int ALIGNMENT>
			struct print_type 
			{
			  static const std::string value(std::string const & mem_accessor = "")
			  {
				return mem_accessor + " " + print_aligned_type<T,ALIGNMENT>::value();
			  }
			};

			template<typename T, unsigned int ALIGNMENT>
			struct print_type<T*, ALIGNMENT> 
			{
			  static const std::string value(std::string const & mem_accessor = "")
			  {
				return mem_accessor + " " +  print_type<T,ALIGNMENT>::value() + "*" ;
			  }
			};
    }
}
#endif
