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
#include "viennacl/tools/shared_ptr.hpp"

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

            template<class T>
            struct SharedPtr2Raw {
                T* operator ()(viennacl::tools::shared_ptr<T> const & value ) const { return value.get(); }
            };


            template<class Base,class Target>
            struct UnsafeBase2Target { Target* operator ()( Base* value ) const { return static_cast<Target*>(value); } };

            void replace_all_occurences(std::string& str, const std::string& oldstr, const std::string& newstr)
            {
              size_t pos = 0;
              while((pos = str.find(oldstr, pos)) != std::string::npos)
              {
                 str.replace(pos, oldstr.length(), newstr);
                 pos += newstr.length();
              }
            }

            template<class T>
            struct is_primitive_type{ enum { value = 0 }; };

            template<> struct is_primitive_type<char>{ enum { value = 1}; };
            template<> struct is_primitive_type<unsigned char>{ enum { value = 1}; };
            template<> struct is_primitive_type<int>{ enum { value = 1}; };
            template<> struct is_primitive_type<unsigned int>{ enum { value = 1}; };
            template<> struct is_primitive_type<short>{ enum { value = 1}; };
            template<> struct is_primitive_type<unsigned short>{ enum { value = 1}; };
            template<> struct is_primitive_type<long>{ enum { value = 1}; };
            template<> struct is_primitive_type<unsigned long>{ enum { value = 1}; };
            template<> struct is_primitive_type<float>{ enum { value = 1}; };
            template<> struct is_primitive_type<double>{ enum { value = 1}; };

		    template<class T>
            struct print_type;

            template<>
            struct print_type<char>
            {
              static const std::string value() { return "char"; }
            };

            template<>
            struct print_type<unsigned char>
            {
              static const std::string value() { return "unsigned char"; }
            };


			template<>
            struct print_type<int>
			{
			  static const std::string value() { return "int"; }
			};

			template<>
            struct print_type<unsigned int>
			{
			  static const std::string value() { return "unsigned int"; }
			};

            template<>
            struct print_type<short>
            {
              static const std::string value() { return "short"; }
            };

            template<>
            struct print_type<unsigned short>
            {
              static const std::string value() { return "unsigned short"; }
            };

			template<>
            struct print_type<long>
			{
			  static const std::string value() { return "long"; }
			};

			template<>
            struct print_type<unsigned long>
			{
              static const std::string value() { return "unsigned long"; }
			};

			template<>
            struct print_type<float>
			{
			  static const std::string value() { return "float"; }
			};

			template<>
            struct print_type<double>
			{
			  static const std::string value() { return "double"; }
			};


    }
}
#endif
