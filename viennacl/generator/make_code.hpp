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
#include <list>

namespace viennacl
{
  namespace generator
  {

    /** @brief Inline code for an expression from scalars
        @tparam T Tree specifying the expression
    */
    template <class T>
    struct make_expression_code
    {
      static std::string value(std::string const & loop_accessor)
      {
        return T::name();
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


    static void replace_all_string(std::string & str, std::string const & from, std::string const & to){
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                 str.replace(start_pos, from.length(), to);
                 start_pos += to.length();
        }
    }

    /**
     * @brief The infos_base class
     */
    class infos_base{
    public:
        infos_base(std::string const & scalartype,
                   std::string const & name): scalartype_(scalartype), name_(name){ }

        std::string const & name() const{ return name_; }
        std::string const & scalartype() const{ return scalartype_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string access_name() const { return access_name_; }
    private:
        std::string scalartype_;
        std::string name_;
        mutable std::string access_name_;
    };

    /**
     * The expr_infos_base class
     */
    template<class T>
    class expr_infos_base{
    public:
        expr_infos_base(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

        void add_infos(T const & infos){
            data_.push_back(infos);
        }

        std::list<T> const & data() const{ return data_; }

        std::string generate() const{
            std::string res(expression_string_base_);
            for(typename std::list<T>::const_iterator it = data_.begin() ; it!= data_.end() ; ++it){
                replace_all_string(res,it->name(),it->access_name());
            }
            return res;
        }

    private:
        std::list<T> data_;
        mutable std::string expression_string_base_;
    };


    /**
     * @brief The mat_infos class
     */
    class mat_infos : public infos_base{
    public:
        mat_infos(std::string const & scalartype,
                  std::string const & name,
                  std::string const & size1,
                  std::string const & size2,
                  bool const & is_rowmajor,
                  bool const & is_transposed) :infos_base(scalartype,name)
                                               ,  size1_(size1)
                                              , size2_(size2)
                                              , is_rowmajor_(is_rowmajor)
                                              , is_transposed_(is_transposed){ }
        std::string const & size1() const{ return size1_; }
        std::string const & size2() const{ return size2_; }
        bool const is_rowmajor() const { return is_rowmajor_; }
        bool const is_transposed() const { return is_transposed_; }
    private:
        std::string size1_;
        std::string size2_;
        bool is_rowmajor_;
        bool is_transposed_;
    };


    /**
     * @brief The vec_infos class
     */
    class vec_infos : public infos_base{
    public:
        vec_infos(std::string const & scalartype,
                          std::string const & name,
                  std::string const & size) : infos_base(scalartype,name), size_(size){ }
        std::string const & size() const{ return size_; }

    private:
        std::string size_;
    };

    /**
     * @brief The scal_infos class
     */
    class scal_infos : public infos_base{
    public:
        scal_infos(std::string const & scalartype, std::string const & name) : infos_base(scalartype,name) { }
    };

    /**
     * @brief The mat_expr_infos class
     */
    typedef expr_infos_base<mat_infos> mat_expr_infos;


    /**
     * @brief The vec_expr_infos class
     */
    typedef expr_infos_base<vec_infos> vec_expr_infos;

    /**
     * @brief The scal_expr_infos class
     */
    typedef expr_infos_base<scal_infos> scal_expr_infos;


    template<class AssignedT, class LhsT, class RhsT, class AddExprT>
    class prod_infos_base : infos_base{
    public:
        prod_infos_base(std::string const & name
                                  ,AssignedT const & assigned
                                  ,std::string const & op_str
                                  ,LhsT const & lhs
                                  ,RhsT const & rhs
                                  ,AddExprT const & additional_expression
                                  ,std::string const & expression_string_base) : infos_base(assigned.scalartype(),name)
                                                                              , assigned_(assigned)
                                                                              , lhs_(lhs)
                                                                              , op_str_(op_str)
                                                                              , rhs_(rhs)
                                                                              , additional_expression_(additional_expression)
                                                                              , expression_string_base_(expression_string_base){ }
        LhsT const & lhs() const { return lhs_; }
        RhsT const & rhs() const { return rhs_; }
        AddExprT const & additional_expression() const { return additional_expression_; }
        AssignedT const & assigned() const { return assigned_; }
        std::string generate(){
            std::string tmp(expression_string_base_);
            replace_all_string(tmp,name_,access_name_);
            return assigned_.access_name() + op_str_ + tmp + "+" + additional_expression_.generate();
        }

    private:
        std::string op_str_;
        AssignedT assigned_;
        LhsT lhs_;
        RhsT rhs_;
        AddExprT additional_expression_;
        mutable std::string expression_string_base_;
    };

    /**
     * @brief The matvec_prod_infos struct
     */
    class matvec_prod_infos : public prod_infos_base<vec_infos, mat_expr_infos, vec_expr_infos, vec_expr_infos>{
        typedef prod_infos_base<vec_infos, mat_expr_infos, vec_expr_infos, vec_expr_infos> BaseT;
    public:
        matvec_prod_infos(std::string const & name
                          ,vec_infos const & assigned
                          ,std::string const & op_str
                          ,mat_expr_infos const & lhs
                          ,vec_expr_infos const & rhs
                          ,vec_expr_infos const & additional_expression
                          ,std::string const & expression_string_base): BaseT(name,assigned,op_str,lhs,rhs,additional_expression,expression_string_base){ }
    };

    /**
     * @brief The inprod_infos struct
     */




    template<class U>
    static vec_infos wrap_vec(){
        return vec_infos(print_type<typename U::ScalarType,1>::value(),
                                 U::name(),
                                 U::internal_size2_name());
    }

    template<class U>
    static mat_infos wrap_mat(){
       return mat_infos(print_type<typename U::ScalarType,1>::value(),
                                      U::name(),
                                      U::size1_name(),
                                      U::size2_name(),
                                      result_of::is_row_major<U>::value,
                                      false);
    }

    template<class U>
    static scal_infos wrap_scal(){
       return scal_infos(print_type<typename U::ScalarType,1>::value(),
                                      U::name());
    }

    template<class U>
    struct foreach_functor{
        static void execute(mat_expr_infos & res){ res.add_infos(wrap_mat<U>()); }
        static void execute(vec_expr_infos & res){ res.add_infos(wrap_vec<U>()); }
        static void execute(scal_expr_infos & res){ res.add_infos(wrap_scal<U>()); }
    };

    /**
     * @brief wrap_vec_expr class
     */
    template<class U>
    static vec_expr_infos wrap_vecexpr(){
            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_vector>::Result List;
            vec_expr_infos res(make_expression_code<U>::value(""));
            typelist_utils::ForEach<List,foreach_functor>::execute(res);
            return res;
    }

    template<class U>
    static mat_expr_infos wrap_matexpr(){
            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_matrix>::Result List;
            mat_expr_infos res(make_expression_code<U>::value(""));
            typelist_utils::ForEach<List,foreach_functor>::execute(res);
            return res;
    }

    template<class U>
    static scal_expr_infos wrap_scalexpr(){
            typedef typename tree_utils::extract_if<U,result_of::is_symbolic_scalar>::Result List;
            scal_expr_infos res(make_expression_code<U>::value(""));
            typelist_utils::ForEach<List,foreach_functor>::execute(res);
            return res;
    }

    template<unsigned int>
    class inprod_infos;

    template<>
    class inprod_infos<1>: public prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos>{
        typedef prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos> BaseT;
    public:
        inprod_infos(std::string const & name
                          ,vec_expr_infos const & lhs
                     ,vec_expr_infos const & rhs): BaseT(name,scal_infos("",""),"",lhs,rhs,wrap_scalexpr<NullType>(),""){ }
    };


    template<>
    class inprod_infos<2>: public prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos>{
        typedef prod_infos_base<scal_infos,vec_expr_infos,vec_expr_infos, scal_expr_infos> BaseT;
    public:
        inprod_infos(std::string const & name
                          ,scal_infos const & assigned
                          ,std::string const & op_str
                          ,vec_expr_infos const & lhs
                          ,vec_expr_infos const & rhs
                          ,scal_expr_infos const & additional_expression
                     ,std::string const & expression_string_base): BaseT(name,assigned,op_str,lhs,rhs,additional_expression,expression_string_base){ }
    };

    template<class T>
    static matvec_prod_infos wrap_matvec_prod(){
        typedef typename tree_utils::extract_if<T,result_of::is_product_leaf>::Result::Head Product;
        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vectors;
        typedef typename Product::LHS ProdLHS;
        typedef typename Product::RHS ProdRHS;
        typedef typename result_of::expression_type<T>::Result    IntermediateType;
        static const unsigned int Alignment = IntermediateType::Alignment;
            VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);
            return matvec_prod_infos(Product::name()
                                     ,wrap_vec<typename T::LHS>()
                                     ,T::OP::expression_string()
                                     ,wrap_matexpr<ProdLHS>()
                                     ,wrap_vecexpr<ProdRHS>()
                                     ,wrap_vecexpr<Vectors>()
                                     ,make_expression_code<Product>::value(""));
    }

    template<class T>
    struct wrap_inprod{

            template<class U>
            struct functor_compute{
                typedef typename U::Type InProd;
                typedef typename InProd::LHS InProdLHS;
                typedef typename InProd::RHS InProdRHS;

                static void execute(std::vector<inprod_infos<1> > & arg){
                    arg.push_back(inprod_infos<1>(InProd::name()
                                                  ,wrap_vecexpr<InProdLHS>()
                                                  ,wrap_vecexpr<InProdRHS>()));
                }
            };

            template<class U>
            struct functor_reduce{
                typedef typename tree_utils::extract_if<T,result_of::is_inner_product_leaf >::Result::Head InProd;
                typedef typename tree_utils::remove_if<T,result_of::is_inner_product_leaf,false >::Result Scalars;
                typedef typename InProd::LHS InProdLHS;
                typedef typename InProd::RHS InProdRHS;
                typedef typename result_of::expression_type<T>::Result    IntermediateType;


                static void execute(std::vector<inprod_infos<2> > & arg){
                    arg.push_back(inprod_infos<2>(InProd::name()
                                                  ,wrap_scal<typename T::LHS>()
                                                  ,T::OP::expression_string()
                                                  ,wrap_vecexpr<InProdLHS>()
                                                  ,wrap_vecexpr<InProdRHS>()
                                                  ,wrap_scalexpr<Scalars>()
                                                  ,make_expression_code<InProd>::value("")));
                }
            };

            static void execute(std::vector<inprod_infos<1> > & arg_compute, std::vector<inprod_infos<2> > & arg_reduce){
                typelist_utils::ForEach< typename tree_utils::extract_if<T,result_of::is_inner_product_impl>::Result, functor_compute>::execute(arg_compute);
                typelist_utils::ForEach< typename tree_utils::extract_if<T,result_of::is_inner_product_leaf>::Result, functor_reduce>::execute(arg_reduce);
            }
    };

    static std::string generate_matvec_prod_backend(std::vector<matvec_prod_infos> const & infos){
        std::string res;
        mat_infos const & mat = infos.front().lhs().data().front();
        vec_infos const & vec = infos.front().rhs().data().front();
        res += "   __local shared_memory_ptr[64];\n";
        res += "   unsigned int row_gid = get_global_id(0)/get_local_size(0);\n";
        res += "   unsigned int col_gid = get_global_id(0)%get_local_size(0);\n";
        res += "   unsigned int lid = get_local_id(0);\n";
        res += "   for(unsigned int row = row_gid ; row < " + mat.size1() + " ; row+=get_num_groups(0)){\n";
        res += "       " + mat.scalartype() + " sum = 0;\n";
        res += "       for(unsigned int col = col_gid ; col < " + mat.size2() + " ; col+=get_local_size(0)){\n";
        mat.access_name(mat.name()+"row");
        vec.access_name("col");
        res += "            sum +=  " + infos.front().lhs().generate() + "*" +  infos.front().rhs().generate() + ";\n";
        res += "       }\n";
        res += "       shared_memory_ptr[lid]=sum;\n";
        res += "       for(unsigned int stride=get_local_size(0)/2 ; stride>0 ; stride>>=1){\n";
        res += "           barrier(CLK_LOCAL_MEM_FENCE);\n";
        res += "           if(lid < stride) shared_memory_ptr[lid]+=shared_memory_ptr[lid+stride];\n";
        res += "       }\n";
        res += "       if(lid==0) shared_memory_ptr[0] ;\n";
        res+=";\n";
        res += "   }\n";
        res += "}\n";
        return res;
    }

    static std::string generate_blas1_backend(std::vector<inprod_infos<2> > & inprods_to_reduce,
                        std::vector<vec_expr_infos> & vec_exprs,
                        std::vector<inprod_infos<1> > & inprods_to_compute){
        return "blabla";

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
                      "   __local shared_memory_ptr[64];\n"
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




    template<class T>
    struct make_code<MatVecToken<T> >
    {

      template<class U>
      struct fill_infos{
          static void execute(std::vector<matvec_prod_infos> & res){
              res.push_back(wrap_matvec_prod<U>::create());
          }
      };

//    /** @brief Functor to make code from a token
//        @tparam TOKEN The token to make the code for
//    */
//    template <class TOKEN>
//    struct make_code;

//    template<>
//    struct make_code<NullType>
//    {
//      static std::string value() { return ""; }

//      static std::string sum() { return ""; }

//      static std::string reduction() { return ""; }
//    };

//    template <class EXPR>
//    struct make_code<ArithmeticToken<EXPR> >
//    {
//      static std::string value()
//      {
//        std::string res;
//        res+="\n//Arithmetic Token\n";
//        res+=make_expression_code<EXPR>::value("gid") + ";\n";
//        return res;
//      }
//    };

//    template <class T>
//    struct make_code<InProdToken<T, 1>  >
//    {
//      template<class U>
//      struct generate_code_sum
//      {
//        private:
//          typedef typename U::Type ARG;
//          typedef typename ARG::LHS LHS;
//          typedef typename ARG::RHS RHS;
//        public:
//          static void execute(std::string & generated_code)
//          {
////              generated_code+= U::private_value() + " += " + dot_product<LHS,RHS>::value("gid","gid") + ";\n";
//          }
//      };


//      template<class U>
//      struct generate_code_reduction
//      {
//        private:
//          typedef typename U::Type ARG;
//          typedef typename ARG::LHS LHS;
//          typedef typename ARG::RHS RHS;
//        public :

//          static void execute(std::string & generated_code)
//          {
//              generated_code+=
//                      "shared_memory_ptr[get_local_id(0)] = " + U::private_value() + ";\n"
//                      "for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
//                      "{\n"
//                      "  barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "  if (get_local_id(0) < stride)\n"
//                      "  shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
//                      "}\n"
//                      +ARG::name() + "[get_group_id(0)] = shared_memory_ptr[0];";
//          }
//      };

//      static std::string sum()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code_sum>::execute(res);
//          return res;
//      }

//      static std::string reduction()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code_reduction>::execute(res);
//          return res;
//      }

//    };

//    template <class T>
//    struct make_code<InProdToken<T, 0> >
//    {
//      template<class U>
//      struct generate_code
//      {
//        private:
//          typedef typename U::LHS LHS;
//          typedef typename U::RHS RHS;

//        public:

//          static void execute(std::string & generated_code)
//          {
//              generated_code+=
//                      "{\n"
//                      "   __local shared_memory_ptr[64];\n"
//                      "   float sum = 0;\n"
//                      "   for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n"
//                      "   {\n"
//                      "      sum+= " +U::name() +"[i];\n"
//                      "   };\n"
//                      "   shared_memory_ptr[get_local_id(0)]=sum;\n"
//                      "   for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
//                      "   {\n"
//                      "      barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "      if (get_local_id(0) < stride)\n"
//                      "      shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
//                      "   }\n"
//                      "   if(get_local_id(0)==0)\n"
//                      "       "+U::local_value() + " = shared_memory_ptr[0];\n"
//                      "   barrier(CLK_LOCAL_MEM_FENCE);\n"
//                      "}\n";
//          }
//      };

//      static std::string value()
//      {
//          std::string res;
//          typelist_utils::ForEach<T,generate_code>::execute(res);
//          return res;
//      }
//    };





//    template<class T>
//    struct make_code<MatVecToken<T> >
//    {

//      template<class U>
//      struct fill_infos{
//          static void execute(std::vector<matvec_prod_infos> & res){
//              res.push_back(wrap_matvec_prod<U>());
//          }
//      };

//      public:
//        static std::string value()
//        {
//          std::vector<matvec_prod_infos> infos;
//          typelist_utils::ForEach<T,fill_infos>::execute(infos);
//          return generate_matvec_prod_backend(infos);
//        }

//    };


  }

}

#endif

