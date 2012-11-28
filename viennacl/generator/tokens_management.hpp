#ifndef VIENNACL_GENERATOR_TOKENS_MANAGEMENT_HPP
#define VIENNACL_GENERATOR_TOKENS_MANAGEMENT_HPP

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

/** @file viennacl/generator/tokens_management.hpp
 *  @brief Creation and management of the tokens list
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/symbolic_types.hpp"
#include "viennacl/generator/operators.hpp"
#include "viennacl/generator/tree_operations.hpp"
#include "viennacl/generator/meta_tools/typelist.hpp"
#include "viennacl/generator/make_code.hpp"
#include "viennacl/generator/tweaking.hpp"

namespace viennacl
{
  namespace generator
  {

    //Tokens

    /** @brief Base structure for representing Token */
    template <class Expr_>
    struct Token
    {
      typedef Expr_ Expr;
    };

    /** @brief Base structure for representing Matrix-Matrix Product token */
    template <class Expr>
    struct MatMatToken : public Token<Expr>
    {

    };

    /** @brief Base structure for representing Matrix-Vector Product token */
    template <class Expr>
    struct MatVecToken : public Token<Expr>
    {
    };

    /** @brief Base structure for representing Inner Product token
        @tparam Step The inner product is calculated in 2 steps with different code/token for each step
    */
    template <class Expr, unsigned int Step_>
    struct InProdToken : public Token<Expr>
    {
        static const int Step = Step_;
    };

    template<class Expr>
    struct ArithmeticToken : public Token<Expr>
    {

    };

//    // Traits

//    namespace result_of
//    {

//      template<class T>
//      struct is_vector_assignment
//      {
//        enum { value = 0 };
//      };

//      template<class LHS, class OP, class RHS>
//      struct is_vector_assignment<compound_node<LHS,OP,RHS> >
//      {
//        enum { value = result_of::is_assignment<OP>::value && result_of::is_symbolic_vector<LHS>::value };
//      };

//      template<class Bound, class Expr>
//      struct is_vector_assignment<repeater_impl<Bound,Expr> >
//      {
//        enum { value = 1 };
//      };

//      template<class T>
//      struct is_scalar_assignment
//      {
//        enum { value = 0 };
//      };

//      template<class LHS, class OP, class RHS>
//      struct is_scalar_assignment<compound_node<LHS,OP,RHS> >
//      {
//        enum { value = result_of::is_assignment<OP>::value && result_of::is_symbolic_gpu_scalar<LHS>::value };
//      };

//    }

//    template<class T, class Enable = void>
//    struct get_operations_lhs
//    {
//        typedef NullType Result;
//    };

//    template<class T>
//    struct get_operations_lhs<T,typename viennacl::enable_if<result_of::is_assignment_compound<T>::value>::type>
//    {
//        typedef typename T::LHS Result;
//    };

//    template<class Head, class Tail>
//    struct get_operations_lhs<typelist<Head,Tail> >
//    {
//        typedef typelist<typename get_operations_lhs<Head>::Result,
//                         typename get_operations_lhs<Tail>::Result> Result;
//    };




    template<class T>
    class fill_blas1_infos{

        template<class U>
        static void exec_impl(std::list<vec_expr_infos_base *> & vec_exprs, std::list<scal_expr_infos_base *> & scal_exprs, typename viennacl::enable_if<result_of::is_vector_expression<U>::value>::type* dummy = 0){
            vec_exprs.push_back(& vec_expr_infos<U>::get());
        }

        template<class U>
        static void exec_impl(std::list<vec_expr_infos_base *> & vec_exprs, std::list<scal_expr_infos_base *> & scal_exprs, typename viennacl::enable_if<result_of::is_scalar_expression<U>::value>::type* dummy = 0){
            scal_exprs.push_back(& scal_expr_infos<U>::get());
        }

    public:

        static void execute(std::list<vec_expr_infos_base *> & vec_exprs, std::list<scal_expr_infos_base *> & scal_exprs){
            exec_impl<T>(vec_exprs,scal_exprs);
        }
    };

    template<class ExpressionsList>
    class body_code
    {
    private:
//        template<class U>
//        struct fill_matvec_prod_infos{
//            static void execute(std::vector<matvec_prod_infos> & res){
//                res.push_back(wrap_matvec_prod<U>());
//            }
//        };

//        template<class U>
//        struct fill_blas1_infos{



//            static void execute(std::vector<inprod_infos<1> > & inprods_compute,std::vector<inprod_infos<2> > & inprods_reduce){
//                wrap_inprod<U>::execute(inprods_compute,inprods_reduce);
//            }

//            static void execute(std::vector<vec_expr_infos> & vec_exprs){
//                bool has_ip1 = tree_utils::count_if<U,result_of::is_inner_product_leaf>::value;
//                bool has_ip2 = tree_utils::count_if<U,result_of::is_inner_product_impl>::value;
//                if(!has_ip1 && !has_ip2){
//                    vec_exprs.push_back(wrap_vecexpr<U>());
//                }
//            }

//        };

//        template<class U>
//        static const std::string value_impl(typename viennacl::enable_if< tree_utils::count_if<U,result_of::is_product_leaf>::value >::type * dummy = 0){
//            std::vector<matvec_prod_infos> infos;
//            typelist_utils::ForEach<ExpressionsList,fill_matvec_prod_infos>::execute(infos);
//            return generate_matvec_prod_backend(infos);
//        }

        template<class U>
        static const std::string value_impl(typename viennacl::enable_if<! tree_utils::count_if<U,result_of::is_product_leaf>::value >::type * dummy = 0){
            std::list<vec_expr_infos_base *> vec_exprs;
            std::list<scal_expr_infos_base *> scal_exprs;
            typelist_utils::ForEach<ExpressionsList,fill_blas1_infos>::execute(vec_exprs,scal_exprs);
            return blas1_generator(vec_exprs,scal_exprs)();
        }

    public:
        /** @brief generates the actual code */
        static const std::string value()
        {
            return value_impl<ExpressionsList>();
        }
    };


//    /** @brief Functor to generates the body code of a kernel from a typelist of expressions.
//        @tparam ExpressionsList Typelist of expressions.
//    */
//    template<class ExpressionsList>
//    struct body_code
//    {
//      private:

//        template<class T>
//        struct requires_vector_access
//        {
//          enum{ value = result_of::is_vector_assignment<T>::value || result_of::is_inner_product_impl<T>::value };
//        };

//        typedef typename get_operations_from_expressions<ExpressionsList>::Result OperationsList;


//        template<class TList, template<class>class Pred>
//        struct fill_expression_updates
//        {

//          template<class Tree>
//          struct functor
//          {
//            private:

//              template<class T>
//              static void execute_impl(T, std::string & generated_code, unsigned int & /*nested_repeats_counter*/)
//              {
//                  if(Pred<Tree>::value)
//                      generated_code+=make_expression_code<Tree>::value("gid") + ";\n";
//              }

//              template<class Bound, class Expressions>
//              static void execute_impl(repeater_impl<Bound,Expressions> ,std::string & generated_code, unsigned int & nested_repeats_counter)
//              {
//                  std::string repeater_name = "Repeat"+to_string(nested_repeats_counter) ;
//                  if(tree_utils::count_if<Expressions,Pred>::value){
//                      generated_code += "for(int " + repeater_name + " = 0 ; " + repeater_name + " < " + Bound::name() + " ; ++" + repeater_name + "){\n";
//                      ++nested_repeats_counter;
//                      typelist_utils::ForEach<Expressions, fill_expression_updates::functor>::execute(generated_code, nested_repeats_counter);
//                      generated_code += "}\n";
//                  }
//              }
              
//            public:
//              static void execute(std::string & generated_code, unsigned int & nested_repeats_counter)
//              {
//                execute_impl(Tree(),generated_code, nested_repeats_counter);
//              }

//          };

//          static void execute(std::string & generated_code)
//          {
//            unsigned int dummy=0;
//            typelist_utils::ForEach<TList,functor>::execute(generated_code,dummy);
//          }
//        };

//        /** @brief Functor to store the values in registers in case of multiple uses
//            @tparam Pred predicate to specify the type of the value to store ( result_of::is_symbolic_scalar, result_of::is_symbolic_vector ...)
//        */
//        template<template<class> class Pred>
//        struct declarations
//        {
//          private:
//            template<class T>
//            struct functor
//            {
//              static void execute(std::string &generated_code)
//              {
//                generated_code+=T::declarations();
//              }
//            };
            
//          public:
//            static void execute(std::string & generated_code)
//            {
//              typedef typename tree_utils::extract_if<OperationsList,Pred>::Result PredValidTmp;
//              typedef typename typelist_utils::no_duplicates<PredValidTmp>::Result PredValid;
//              typelist_utils::ForEach<PredValid,functor>::execute(generated_code);
//            }
//        };


//        /** @brief Functor to store the values back in registers in case of multiple uses
//            @tparam Pred predicate to specify the type of the value to store ( result_of::is_symbolic_scalar, result_of::is_symbolic_vector ...)
//        */
//        template<template<class> class Pred>
//        struct assignements
//        {
//          private:
//            template<class T>
//            struct functor
//            {
//              static void execute(std::string &generated_code)
//              {
//                generated_code += T::assignements();
//              }
//            };
            
//          public:
//            static void execute(std::string & generated_code)
//            {
//              typelist_utils::ForEach<typename typelist_utils::no_duplicates< typename tree_utils::extract_if<typename  get_operations_lhs<OperationsList>::Result
//                                                                                                              ,Pred >::Result
//                                                                             >::Result,
//                                      functor>::execute(generated_code);
//            }
//        };

//        /** @brief Generates code for vector expressions, if any
//            @param cond Condition for the static_if
//            @param res Reference to the result
//        */
//        template<class RequireGidLoop>
//        static void fill_vector_expression(Int2Type<false> /*cond*/, std::string & res)
//        {
//          typedef typename get_operations_from_expressions<RequireGidLoop>::Result RequireGidOperations;
//          typedef typename tree_utils::extract_if<RequireGidOperations,result_of::is_symbolic_vector>::Result SymVecs;
//          typedef typename result_of::expression_type<typename SymVecs::Head>::Result VecExpr;
//          std::string bound = VecExpr::internal_size_expression();
//          res += "for(unsigned int gid=get_global_id(0) ; gid < " + bound + " ; gid+=get_global_size(0))\n";
//          res += "{\n";
//          //For each unique symbolic vector or symbolic matrix in the tree, store the gid value in a local register
//          declarations<result_of::or_is<result_of::is_symbolic_vector,result_of::is_symbolic_matrix>::Pred>::execute(res);
//          res += "\n";

//          fill_expression_updates<RequireGidLoop,result_of::is_vector_assignment>::execute(res);

//          //Inner Product - Step 1 - Sum
//          typedef typename tree_utils::extract_if_unique<RequireGidLoop,result_of::is_inner_product_impl>::Result InProd;
//          typedef typename get_type_if<NullType,InProdToken<InProd,1>,result_of::is_null_type<InProd>::value>::Result InProdToken_;
//          res += make_code<InProdToken_>::sum();

//          assignements<result_of::is_symbolic_vector>::execute(res);
//          res += "}\n";
//        }

//        /** @brief Does not generate anything */
//        template<class RequireGidLoop>
//        static void fill_vector_expression(Int2Type<true>, std::string &) {}

//        /** @brief Generates code for simple vector expressions */
//        static const std::string vector_code_impl(Int2Type<0> /*linear_expression*/)
//        {
//          typedef typename tree_utils::extract_if<ExpressionsList,requires_vector_access>::Result RequireGidLoopTmp;
//          typedef typename typelist_utils::no_duplicates<RequireGidLoopTmp>::Result RequireGidLoop;
//          std::string res;
//          fill_vector_expression<RequireGidLoop>(Int2Type<result_of::is_null_type<RequireGidLoop>::value>(),res);
//          return res;
//        }

//        /** @brief Generates code for simple matrix-vector product */
//        static const std::string vector_code_impl(Int2Type<1> /*matvec_prod*/)
//        {
//          std::string res;
//          res += make_code<MatVecToken<OperationsList> >::value();
//          return res;
//        }

//      public:

//        /** @brief generates the actual code */
//        static const std::string value()
//        {
//          std::string res;
//          res += "{\n";


//          declarations<result_of::is_symbolic_gpu_scalar>::execute(res);
//          declarations<result_of::is_inner_product_impl>::execute(res);
//          declarations<result_of::is_inner_product_leaf>::execute(res);

//          //Inner Product - Step 2 - Final Reduction
//          typedef typename tree_utils::extract_if_unique<OperationsList,result_of::is_inner_product_leaf>::Result InProd0;
//          typedef typename get_type_if<NullType,InProdToken<InProd0,0>,result_of::is_null_type<InProd0>::value>::Result InProd0Token_;
//          res+=make_code<InProd0Token_>::value();

//          if(tree_utils::count_if<OperationsList, requires_vector_access>::value)
//              res += vector_code_impl(Int2Type< (tree_utils::count_if<OperationsList, result_of::is_product_leaf>::value > 0) >());

//          //Inner Product - Step 1 - Reduction
//          typedef typename tree_utils::extract_if_unique<OperationsList,result_of::is_inner_product_impl>::Result InProd1;
//          typedef typename get_type_if<NullType,InProdToken<InProd1,1>,result_of::is_null_type<InProd1>::value>::Result InProd1Token_;
//          res += make_code<InProd1Token_>::reduction();

//          fill_expression_updates<ExpressionsList,result_of::is_scalar_assignment>::execute(res);
//          assignements<result_of::is_symbolic_gpu_scalar>::execute(res);

//          res += "}\n";
//          return res;
//        }
//    };
  }
}
#endif
