#ifndef VIENNACL_GENERATOR_MAKE_CODE_HPP
#define VIENNACL_GENERATOR_MAKE_CODE_HPP

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


/** @file viennacl/generator/make_code.hpp
 *  @brief Definition of code generation policies. Experimental.
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/forwards.h"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/symbolic_types.hpp"
#include "viennacl/generator/result_of.hpp"
#include "viennacl/generator/tree_operations.hpp"

namespace viennacl
{
  namespace generator
  {



    template <class T>
    struct inner_prod_impl_t
    {
      typedef T Type;

      static std::string name()
      {
        return T::name();
      }

      static std::string private_value()
      {
        return "private_"+name();
      }

      static std::string declarations()
      {
        return print_type<typename T::ScalarType,1>::value() + " " + private_value() + "=0;\n" ;
      }

      static std::string kernel_arguments()
      {
        return T::kernel_arguments();
      }

      enum { id = T::id };
    };

    /** @brief Inline code for an expression from scalars
        @tparam T Tree specifying the expression
    */
    template <class T>
    struct make_expression_code
    {
      static std::string value(std::string const & loop_accessor)
      {
        if(loop_accessor == "gid") return  T::gid_val_name();
        else return T::name() + '[' + loop_accessor + ']';
      }
    };

    template <long VAL>
    struct make_expression_code<symbolic_constant<VAL> >
    {
      static std::string value(std::string const & )
      {
        return to_string(VAL);
      }
    };

    template<class T>
    struct make_expression_code<inner_prod_impl_t<T> >
    {
      static std::string value(std::string const & )
      {
        return "";
      }
    };

    template <unsigned int ID, class SCALARTYPE>
    struct make_expression_code<cpu_symbolic_scalar<ID,SCALARTYPE> >
    {
      static std::string value(std::string const & )
      {
        return  cpu_symbolic_scalar<ID,SCALARTYPE>::name();
      }
    };

    template <unsigned int ID,class SCALARTYPE>
    struct make_expression_code<gpu_symbolic_scalar<ID,SCALARTYPE> >
    {
      static std::string value(std::string const & )
      {
        return  gpu_symbolic_scalar<ID,SCALARTYPE>::val_name();
      }
    };

    template <class LHS, class RHS>
    struct make_expression_code<compound_node<LHS,inner_prod_type,RHS > >
    {
      private:
        typedef compound_node<LHS,inner_prod_type,RHS> T;

      public:
        static std::string value(std::string const & )
        {
          return T::local_value();
        }
    };


    template <class LHS, class RHS>
    struct make_expression_code<compound_node<LHS,prod_type,RHS > >
    {
      private:
        typedef compound_node<LHS,prod_type,RHS> T;

      public:
        static std::string value(std::string const & )
        {
            return "__PRODVAL__";
        }
    };

    template< >
    struct make_expression_code< NullType >
    {
      static std::string value(std::string const & )
      {
        return "";
      }
    };

    template<class T>
    struct make_expression_code< elementwise_modifier<T> >
    {
      static std::string value ( std::string const & loop_accessor)
      {
        return elementwise_modifier<T>::modify(make_expression_code<typename T::PRIOR_TYPE>::value(loop_accessor));
      }
    };

    template<class LHS, class OP, class RHS >
    struct make_expression_code<compound_node<LHS, OP, RHS> >
    {
      static const std::string value(std::string const & loop_accessor)
      {
        return " ( " + make_expression_code<LHS>::value(loop_accessor)
                + OP::expression_string()
                + make_expression_code<RHS>::value(loop_accessor) + " ) ";
      }
    };



    template<class LHS, class RHS, unsigned int Alignment>
    struct dot_product_impl
    {
      static std::string value(std::string lhs_loop_id,
                               std::string rhs_loop_id)
      {
        return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + "," + make_expression_code<RHS>::value(rhs_loop_id) + ")";
      }
    };

    template<class LHS, class RHS>
    struct dot_product_impl<LHS, RHS, 8>
    {
      static std::string value(std::string lhs_loop_id,
                                std::string rhs_loop_id)
      {
          return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + ".s0123" + ","
                  + make_expression_code<RHS>::value(rhs_loop_id) + ".s0123 )"
                  +  " + dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s4567" + ","
                  + make_expression_code<RHS>::value(rhs_loop_id) + ".s4567 );"
                  ;
      }
    };

    template<class LHS, class RHS>
    struct dot_product_impl<LHS, RHS, 16>
    {
      static std::string value(std::string lhs_loop_id,std::string rhs_loop_id)
      {
        return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + ".s0123" + ","
                + make_expression_code<RHS>::value(rhs_loop_id) + ".s0123)"
                +"\n + dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s4567" + ","
                + make_expression_code<RHS>::value(rhs_loop_id) + ".s4567) "
                +"\n + dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s89ab" + ","
                + make_expression_code<RHS>::value ( rhs_loop_id ) + ".s89ab) "
                +"\n + dot("	+ make_expression_code<LHS>::value ( lhs_loop_id ) + ".scdef" + ","
                + make_expression_code<RHS>::value ( rhs_loop_id ) + ".scdef)"
                ;
      }
    };

    template<class LHS, class RHS>
    struct dot_product
    {
      static std::string value(std::string lhs_loop_id,std::string rhs_loop_id)
      {
        return dot_product_impl<LHS,RHS,LHS::Alignment>::value(lhs_loop_id,rhs_loop_id);
      }
    };

    template<class LHS>
    struct dot_product<LHS, symbolic_constant<1> >
    {
      static std::string value(std::string lhs_loop_id,std::string rhs_loop_id)
      {
        return make_expression_code<LHS>::value(lhs_loop_id);
      }
    };


    /** @brief Functor to make code from a token
        @tparam TOKEN The token to make the code for
    */
    template <class TOKEN>
    struct make_code;

    template<>
    struct make_code<NullType>
    {
      static std::string value() { return ""; }

      static std::string sum() { return ""; }

      static std::string reduction() { return ""; }
    };

    template <class EXPR>
    struct make_code<ArithmeticToken<EXPR> >
    {
      static std::string value()
      {
        std::string res;
        res+="\n//Arithmetic Token\n";
        res+=make_expression_code<EXPR>::value("gid") + ";\n";
        return res;
      }
    };

    template <class T>
    struct make_code<InProdToken<T, 1>  >
    {
      template<class U>
      struct generate_code_sum
      {
        private:
          typedef typename U::Type ARG;
          typedef typename ARG::LHS LHS;
          typedef typename ARG::RHS RHS;
        public:
          static void execute(std::string & generated_code)
          {
              generated_code+= U::private_value() + " += " + dot_product<LHS,RHS>::value("gid","gid") + ";\n";
          }
      };


      template<class U>
      struct generate_code_reduction
      {
        private:
          typedef typename U::Type ARG;
          typedef typename ARG::LHS LHS;
          typedef typename ARG::RHS RHS;
        public :

          static void execute(std::string & generated_code)
          {
              generated_code+=
                      "shared_memory_ptr[get_local_id(0)] = " + U::private_value() + ";\n"
                      "for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
                      "{\n"
                      "  barrier(CLK_LOCAL_MEM_FENCE);\n"
                      "  if (get_local_id(0) < stride)\n"
                      "  shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
                      "}\n"
                      +ARG::name() + "[get_group_id(0)] = shared_memory_ptr[0];";
          }
      };

      static std::string sum()
      {
          std::string res;
          typelist_utils::ForEach<T,generate_code_sum>::execute(res);
          return res;
      }

      static std::string reduction()
      {
          std::string res;
          typelist_utils::ForEach<T,generate_code_reduction>::execute(res);
          return res;
      }

    };

    template <class T>
    struct make_code<InProdToken<T, 0> >
    {
      template<class U>
      struct generate_code
      {
        private:
          typedef typename U::LHS LHS;
          typedef typename U::RHS RHS;

        public:

          static void execute(std::string & generated_code)
          {
              generated_code+=
                      "{\n"
                      "   float sum = 0;\n"
                      "   for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n"
                      "   {\n"
                      "      sum+= " +U::name() +"[i];\n"
                      "   };\n"
                      "   shared_memory_ptr[get_local_id(0)]=sum;\n"
                      "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
                      "   {\n"
                      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
                      "      if (get_local_id(0) < stride)\n"
                      "      shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
                      "   }\n"
                      "   if(get_local_id(0)==0)\n"
                      "       "+U::local_value() + " = shared_memory_ptr[0];\n"
                      "   barrier(CLK_LOCAL_MEM_FENCE);\n"
                      "}\n";
          }
      };

      static std::string value()
      {
          std::string res;
          typelist_utils::ForEach<T,generate_code>::execute(res);
          return res;
      }
    };

    void replace_all_string(std::string & str, std::string const & from, std::string const & to){
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                 str.replace(start_pos, from.length(), to);
                 start_pos += to.length();
        }
    }

    template<class T, class OP, class Assigned>
    struct make_code<MatVecToken<T,OP,Assigned> >
    {
      private:

        typedef typename tree_utils::extract_if<T,result_of::is_product_leaf>::Result::Head Products;
        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vectors;

        typedef typename result_of::expression_type<T>::Result    IntermediateType;
        static const unsigned int Alignment = IntermediateType::Alignment;

        std::string const & size1() const{ return size1_; }
        std::string const & size2() const{ return size2_; }
        std::string const & name() const{ return name_; }
        std::string const & scalartype() const{ return scalartype_; }
        bool const is_rowmajor() const { return is_rowmajor_; }
        bool const is_transposed() const { return is_transposed_; }

    private:
        std::string scalartype_;
        std::string name_;
        std::string size1_;
        std::string size2_;
        bool is_rowmajor_;
        bool is_transposed_;
        mutable std::string expression_string_base_;
    };

    class vec_infos{
    public:
        vec_infos(std::string const & scalartype,
                          std::string const & name,
                          std::string const & size,
                    std::string const & expression_string_base) : scalartype_(scalartype), name_(name), size_(size), expression_string_base_(expression_string_base) { }

        std::string access_element(std::string const & offset) const{
            std::string res(expression_string_base_);
            replace_all_string(res,"__OFFSET__",offset);
            return res;
        }

        std::string const & scalartype(){ return scalartype_; }
        std::string const & name(){ return name_; }
        std::string const & size(){ return size_; }

    private:
        std::string scalartype_;
        std::string name_;
        std::string size_;
        mutable std::string expression_string_base_;
    };

    struct matvec_prod_infos{
        matvec_prod_infos(mat_infos const & lhs, vec_infos const & rhs, std::string const & assignment_str): lhs_(lhs)
                                                                                                            , rhs_(rhs)
                                                                                                            ,assignment_str_(assignment_str){ }
        std::string access_element(std::string const & prodval_name, std::string const & offset) const{
            std::string res(assignment_str_);
            replace_all_string(res,"__OFFSET__",offset);
            replace_all_string(res,"__PRODVAL__",prodval_name);
            return res;
        }

        mat_infos const & lhs() const { return lhs_; }
        vec_infos const & rhs() const { return rhs_; }

    private:
        mat_infos lhs_;
        vec_infos rhs_;
        mutable std::string assignment_str_;
    };

    template<class T>
    matvec_prod_infos wrap_matvec_prod(){
        typedef typename tree_utils::remove_if<T,result_of::is_symbolic_vector,false >::Result Product;
        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vector;
        typedef typename Product::LHS ProdLHS;
        typedef typename Product::RHS ProdRHS;
        typedef typename result_of::expression_type<T>::Result    IntermediateType;
        static const unsigned int Alignment = IntermediateType::Alignment;
        VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);

        mat_infos lhs_infos(print_type<typename ProdLHS::ScalarType,1>::value(),
                            ProdLHS::name(),
                            ProdLHS::internal_size1_name(),
                            ProdLHS::internal_size2_name(),
                            result_of::is_row_major<ProdLHS>::value,
                            false,
                            make_expression_code<ProdLHS>::value("__OFFSET__"));
        vec_infos rhs_infos(print_type<typename ProdRHS::ScalarType,1>::value(),
                            ProdRHS::name(),
                            ProdRHS::internal_size2_name(),
                            make_expression_code<ProdRHS>::value("__OFFSET__"));

        return matvec_prod_infos(lhs_infos,rhs_infos,make_expression_code<T>::value("__OFFSET__"));
    }

    std::string generate_matvec_prod_backend(std::vector<matvec_prod_infos> const & infos){
        std::string res;
        if (infos.front().lhs().is_rowmajor()){
          res += "{\n";
          res += "   unsigned int row_gid = get_global_id(0)/get_local_size(0);\n";
          res += "   unsigned int col_gid = get_global_id(0)%get_local_size(0);\n";
          res += "   unsigned int lid = get_local_id(0);\n";
          res += "   for(unsigned int row = row_gid ; row < " + infos.front().lhs().size1() + " ; row+=get_num_groups(0)){\n";
          res += "       " + infos.front().lhs().scalartype() + " sum = 0;\n";
          res += "       for(unsigned int col = col_gid ; col < " + infos.front().lhs().size2() + " ; col+=get_local_size(0)){\n";
          res += "            sum +=  " + infos.front().lhs().access_element("row") + "*" +  infos.front().rhs().access_element("col") + ";\n";
          res += "       }\n";
          res += "       shared_memory_ptr[lid]=sum;\n";
          res += "       for(unsigned int stride=get_local_size(0)/2 ; stride>0 ; stride>>=1){\n";
          res += "           barrier(CLK_LOCAL_MEM_FENCE);\n";
          res += "           if(lid < stride) shared_memory_ptr[lid]+=shared_memory_ptr[lid+stride];\n";
          res += "       }\n";
          res += "       if(lid==0) " + infos.front().access_element( " shared_memory_ptr[0] ", "row") + ";\n";
          res+=";\n";
          res += "   }\n";
          res += "}\n";
        }
//        else
//        {
//          res += "   for(unsigned int row = get_global_id(0) ; row < " + FirstMatrix::internal_size1_name() + " ; row+=get_global_size(0)){\n";
//          res += "       " + print_type<typename FirstMatrix::ScalarType,1>::value() + " sum = 0;\n";
//          res += "       for(unsigned int col = 0 ; col < " + FirstMatrix::internal_size2_name() + " ; col++){\n";
//          res += "            sum +=  " + make_expression_code<Products>::value("") + ";\n";
//          res += "       }\n";
//          res += "  }\n";
//        }
        return res;
    }



    template<class T>
    struct make_code<MatVecToken<T> >
    {

      template<class U>
      struct fill_infos{
          static void execute(std::vector<matvec_prod_infos> & res){
              res.push_back(wrap_matvec_prod<U>());
          }
      };

      public:
        static std::string value()
        {
          std::vector<matvec_prod_infos> infos;
          typelist_utils::ForEach<T,fill_infos>::execute(infos);
          return generate_matvec_prod_backend(infos);
        }

    };


//    //Matrix-Matrix Product.

//    template<class T>
//    struct make_code<MatMatToken<T> >
//    {
//      private:

//        static bool replace(std::string& str, const std::string& from, const std::string& to)
//        {
//          size_t start_pos = str.find(from);
//          if(start_pos == std::string::npos)
//              return false;
//          str.replace(start_pos, from.length(), to);
//          return true;
//        }


//        typedef typename tree_utils::remove_if<T, result_of::is_product_leaf>::Result                           SCALAR_EXPR;
//        typedef typename tree_utils::extract_if<T,result_of::is_product_leaf>::Result::Head                     ARG;
//        typedef typename ARG::LHS                                   LHS;
//        typedef typename ARG::RHS                                   RHS;
//        typedef typename result_of::expression_type<LHS>::Result    MatExpr;
//        typedef typename MatExpr::ScalarType                        ScalarType;
//        typedef typename MatExpr::Layout                            Layout;
//        static const unsigned int Alignment = result_of::expression_type<LHS>::Result::Alignment;

//      public:
//        static std::string value()
//        {
//            VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);
//            static const std::size_t block_size = 16;
//            static const std::size_t vector_size =  4;
//            std::string res;
//            res += "{\n";
//            res+= "  size_t row_block_id = get_group_id(1);\n";    //refers to the row index in op(A), op(B)
//            res += "  size_t col_block_id = get_group_id(0);\n";    //refers to the col index in op(A), op(B)
//            res += "  size_t row_thread_id = get_local_id(1);\n";
//            res += "  size_t col_thread_id = get_local_id(0);\n";
//            res += "  __local float As[" + to_string(block_size * block_size) + "];\n";
//            res += "  float cv[" + to_string(block_size) + "] = {";;
//            for (std::size_t i=0; i<block_size-1; ++i)
//                res += "0,";
//            res += "0};\n" ;

//            if (result_of::is_row_major<LHS>::value && result_of::is_transposed<LHS>::value)
//            {
//                res += "  size_t aBegin = (row_block_id * " + to_string(block_size) + " * " + LHS::col_inc_name() + " + " + LHS::col_start_name() + ") + " + LHS::row_start_name() + " * " + LHS::internal_size2_name()  + ";\n";
//                res += "  size_t aStep = " + to_string(block_size) + " * " + LHS::internal_size2_name()  + " * " + LHS::row_inc_name() + ";\n";
//                res += "  size_t aEnd = aBegin + " + LHS::internal_size2_name()  + " * " + LHS::row_inc_name() + " * " + LHS::size1_name() + ";\n";
//            }
//            else if (result_of::is_row_major<LHS>::value && !result_of::is_transposed<LHS>::value)
//            {
//                res += "  size_t aBegin = (row_block_id * " + to_string(block_size) + " * " + LHS::row_inc_name() + " + " + LHS::row_start_name() + ") * " + LHS::internal_size2_name()  + " + " + LHS::col_start_name() + ";\n";
//                res += "  size_t aStep = " + to_string(block_size) + " * " + LHS::col_inc_name() + ";\n";
//                res += "  size_t aEnd = aBegin + " + LHS::col_inc_name() + " * " + LHS::size2_name() + ";\n";
//            }
//            else if (!result_of::is_row_major<LHS>::value && result_of::is_transposed<LHS>::value)
//            {
//                res += "  size_t aBegin = (row_block_id * " + to_string(block_size) + " * " + LHS::col_inc_name() + " + " + LHS::col_start_name() + ") * " + LHS::internal_size1_name()  + " + " + LHS::row_start_name() + ";\n";
//                res += "  size_t aStep = " + to_string(block_size) + " * " + LHS::row_inc_name() + ";\n";
//                res += "  size_t aEnd = aBegin + " + LHS::row_inc_name() + " * " + LHS::size1_name() + ";\n";
//            }
//            else if (!result_of::is_row_major<LHS>::value && !result_of::is_transposed<LHS>::value)
//            {
//                res += "  size_t aBegin = (row_block_id * " + to_string(block_size) + " * " + LHS::row_inc_name() + " + " + LHS::row_start_name() + ") + " + LHS::col_start_name() + " * " + LHS::internal_size1_name()  + ";\n";
//                res += "  size_t aStep = " + to_string(block_size) + " * " + LHS::internal_size1_name()  + " * " + LHS::col_inc_name() + ";\n";
//                res += "  size_t aEnd = aBegin + " + LHS::internal_size1_name()  + " * " + LHS::col_inc_name() + " * " + LHS::size2_name() + ";\n";
//            }


//            if (result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//            {
//                res += "  size_t bBegin = (col_block_id * " + to_string(block_size * vector_size) + " * " + RHS::row_inc_name() + " + " + RHS::row_start_name() + ") * " + RHS::internal_size2_name()  + " + " + RHS::col_start_name() + ";\n";
//                res += "  size_t bStep = " + to_string(block_size) + " * " + RHS::col_inc_name() + ";\n";
//            }
//            else if (result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//            {
//                res += "  size_t bBegin = (col_block_id * " + to_string(block_size * vector_size) + " * " + RHS::col_inc_name() + " + " + RHS::col_start_name() + ") + " + RHS::row_start_name() + " * " + RHS::internal_size2_name()  + ";\n";
//                res += "  size_t bStep = " + to_string(block_size) + " * " + RHS::row_inc_name() + " * " + RHS::internal_size2_name()  + ";\n";
//            }
//            else if (!result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//            {
//                res += "  size_t bBegin = (col_block_id * " + to_string(block_size * vector_size) + " * " + RHS::row_inc_name() + " + " + RHS::row_start_name() + ") + " + RHS::col_start_name() + " * " + RHS::internal_size1_name()  + ";\n";
//                res += "  size_t bStep = " + to_string(block_size) + " * " + RHS::col_inc_name() + " * " + RHS::internal_size1_name()  + ";\n";
//            }
//            else if (!result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//            {
//                res += "  size_t bBegin = (col_block_id * " + to_string(block_size * vector_size) + " * " + RHS::col_inc_name() + " + " + RHS::col_start_name() + ") * " + RHS::internal_size1_name()  + " + " + RHS::row_start_name() + ";\n";
//                res += "  size_t bStep = " + to_string(block_size) + " * " + RHS::row_inc_name() + ";\n";
//            }

//            res += "  for(size_t a = aBegin, b = bBegin; a < aEnd; a += aStep, b += bStep) { \n";

//            // copy blocks of op(A) to shared memory (op(A) is column-major in shared memory then)
//            res += "    for(size_t i = 0; i < " + to_string(vector_size) + "; i++)  \n";
//            if (result_of::is_row_major<LHS>::value && result_of::is_transposed<LHS>::value)
//                res += "      As[ (i*" + to_string(vector_size) + " + row_thread_id) + " + to_string(block_size) + " * col_thread_id] = (" + make_expression_code<LHS>::value("a + " + LHS::col_inc_name() + " * (i * " + to_string(vector_size) + " + row_thread_id) + " + LHS::internal_size2_name()  + " * " + LHS::row_inc_name() + " * col_thread_id") + ");\n";
//            else if (result_of::is_row_major<LHS>::value && !result_of::is_transposed<LHS>::value)
//                res += "      As[ (i*" + to_string(vector_size) + " + row_thread_id) + " + to_string(block_size) + " * col_thread_id] = (" + make_expression_code<LHS>::value("a + " + LHS::internal_size2_name()  + " * " + LHS::row_inc_name() + " * (i * " + to_string(vector_size) + " + row_thread_id) + " + LHS::col_inc_name() + " * col_thread_id") + ");\n";
//            else if (!result_of::is_row_major<LHS>::value && result_of::is_transposed<LHS>::value)
//                res += "      As[ (i*" + to_string(vector_size) + " + row_thread_id) + " + to_string(block_size) + " * col_thread_id] = (" + make_expression_code<LHS>::value("a + " + LHS::internal_size1_name()  + " * " + LHS::col_inc_name() + " * (i * " + to_string(vector_size) + " + row_thread_id) + " + LHS::row_inc_name() + " * col_thread_id") + ");\n";
//            else if (!result_of::is_row_major<LHS>::value && !result_of::is_transposed<LHS>::value)
//                res += "      As[ (i*" + to_string(vector_size) + " + row_thread_id) + " + to_string(block_size) + " * col_thread_id] = (" + make_expression_code<LHS>::value("a + " + LHS::row_inc_name() + " * (i * " + to_string(vector_size) + " + row_thread_id) + " + LHS::internal_size1_name()  + " * " + LHS::col_inc_name() + " * col_thread_id") + ");\n";
//            res += "\n";
//            res += "    barrier(CLK_LOCAL_MEM_FENCE); \n";

//            // initialize memory pointers
//            res += "\n";
//            res += "    __local  float *ap = As; \n";
//            if (result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//                res += "    __global float *bp = " + RHS::name() + " + (b + (" + to_string(block_size) + " * row_thread_id + col_thread_id) * " + RHS::row_inc_name() + " * " + RHS::internal_size2_name()  + "); \n";
//            else if (result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//                res += "    __global float *bp = " + RHS::name() + " + (b + (" + to_string(block_size) + " * row_thread_id + col_thread_id) * " + RHS::col_inc_name() + "); \n";
//            else if (!result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//                res += "    __global float *bp = " + RHS::name() + " + (b + (" + to_string(block_size) + " * row_thread_id + col_thread_id) * " + RHS::row_inc_name() + "); \n";
//            else if (!result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//                res += "    __global float *bp = " + RHS::name() + " + (b + (" + to_string(block_size) + " * row_thread_id + col_thread_id) * " + RHS::col_inc_name() + " * " + RHS::internal_size1_name()  + "); \n";
//            res += "\n";

//            std::string rhs_expr;

//            if (result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//                rhs_expr = make_expression_code<RHS>::value("i");
//            else if (result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//                rhs_expr = make_expression_code<RHS>::value("i * " + RHS::internal_size2_name());
//            else if (!result_of::is_row_major<RHS>::value && result_of::is_transposed<RHS>::value)
//                rhs_expr = make_expression_code<RHS>::value("i * " + RHS::internal_size1_name());
//            else if (!result_of::is_row_major<RHS>::value && !result_of::is_transposed<RHS>::value)
//                rhs_expr = make_expression_code<RHS>::value("i");



//            replace(rhs_expr,RHS::name(),"bp");

//            // compute
//            res += "    for(size_t i = 0; i < " + to_string(block_size) + "; i++) { \n";
//            res += "      float bv = " + rhs_expr + "; \n";
//            res += "\n";
//            res += "      for(size_t k = 0; k < " + to_string(block_size) + "; k++)  \n";
//            res += "	    cv[k] += ap[k] * bv; \n";
//            res += "\n";
//            res += "      ap += " + to_string(block_size) + "; \n";
//            res += "    } \n";
//            res += "\n";
//            res += "    barrier(CLK_LOCAL_MEM_FENCE); \n";
//            res += "  } \n";

//            // write to C
//            if (result_of::is_row_major<Assigned>::value)
//            {
//                res += "  int c = " + Assigned::internal_size2_name()  + " * (" + Assigned::row_inc_name() + " * " + to_string(block_size) + " * row_block_id + " + Assigned::row_start_name() + ") + "  //block row index
//                        + to_string(vector_size * block_size) + " * " + Assigned::col_inc_name() + " * col_block_id + " + Assigned::col_start_name() + " \n";  //block column index
//                res += "          + " + Assigned::col_inc_name() + " * (" + to_string(block_size) + " * row_thread_id + col_thread_id); \n";
//            }
//            else
//            {
//                res += "  int c = " + Assigned::row_inc_name() + " * " + to_string(block_size) + " * row_block_id + " + Assigned::row_start_name() + " + ("   // block row index
//                        + to_string(vector_size * block_size) + " * " + Assigned::col_inc_name() + " * col_block_id + " + Assigned::col_start_name() + ") * " + Assigned::internal_size1_name()  + " \n";   // block column index
//                res += "          + " + Assigned::internal_size1_name()  +  " * " + Assigned::col_inc_name() + " * (" + to_string(block_size) + " * row_thread_id + col_thread_id); \n";
//            }

//            res += "  for(size_t i = 0; i < " + to_string(block_size) + "; i++) { \n";


//            if (result_of::is_row_major<Assigned>::value)
//            {
//                if(result_of::is_null_type<SCALAR_EXPR>::value){
//                    res += "    " + Assigned::name() + "[c]" + OP::expression_string() + "cv[i]; \n";
//                    res += "      c += " + Assigned::internal_size2_name()  + " * " + Assigned::row_inc_name() + "; \n";
//                }
//                else{
//                    res += "    " + Assigned::name() + "[c]" + OP::expression_string() + make_expression_code<SCALAR_EXPR>::value ("") + "* cv[i]; \n";
//                    res += "      c += " + Assigned::internal_size2_name()  + " * " + Assigned::row_inc_name() + "; \n";
//                }
//            }
//            else
//            {
//                if(result_of::is_null_type<SCALAR_EXPR>::value){
//                    res += "    " + Assigned::name() + "[c]" + OP::expression_string() + "cv[i]; \n";
//                    res += "      c += " + Assigned::row_inc_name() + "; \n";
//                }
//                else{
//                    res += "    " + Assigned::name() + "[c]" + OP::expression_string() + make_expression_code<SCALAR_EXPR>::value ("") + "* cv[i]; \n";
//                    res += "      c += " + Assigned::row_inc_name() + "; \n";
//                }
//            }

//            res += "  } \n";
//            res += "} \n";

//            return res;
//        }
//    };

  }

}

#endif

