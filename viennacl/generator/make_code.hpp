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
     * @brief The mat_infos class
     */
    class mat_infos{
    public:
        mat_infos(std::string const & scalartype,
                  std::string const & name,
                  std::string const & size1,
                  std::string const & size2,
                  bool const & is_rowmajor,
                  bool const & is_transposed) : scalartype_(scalartype), name_(name)
                                                              ,  size1_(size1), size2_(size2)
                                                              , is_rowmajor_(is_rowmajor), is_transposed_(is_transposed){ }
        std::string const & size1() const{ return size1_; }
        std::string const & size2() const{ return size2_; }
        std::string const & name() const{ return name_; }
        std::string const & scalartype() const{ return scalartype_; }
        bool const is_rowmajor() const { return is_rowmajor_; }
        bool const is_transposed() const { return is_transposed_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string access_name() const { return access_name_; }
    private:
        mutable std::string access_name_;
        std::string scalartype_;
        std::string name_;
        std::string size1_;
        std::string size2_;
        bool is_rowmajor_;
        bool is_transposed_;
    };

    /**
     * @brief The mat_expr_infos class
     */
    class mat_expr_infos{
    public:
        mat_expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

        void add_infos(mat_infos const & infos){
            mats_.push_back(infos);
        }

        std::list<mat_infos> const & matrices() const{ return mats_; }

        std::string generate() const{
            std::string res(expression_string_base_);
            for(std::list<mat_infos>::const_iterator it = mats_.begin() ; it!= mats_.end() ; ++it){
                replace_all_string(res,it->name(),it->access_name());
            }
            return res;
        }

    private:
        std::list<mat_infos> mats_;
        mutable std::string expression_string_base_;
    };

    /**
     * @brief The vec_infos class
     */
    class vec_infos{
    public:
        vec_infos(std::string const & scalartype,
                          std::string const & name,
                          std::string const & size) : scalartype_(scalartype), name_(name), size_(size){ }
        std::string const & scalartype() const{ return scalartype_; }
        std::string const & name() const{ return name_; }
        std::string const & size() const{ return size_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string access_name() const { return access_name_; }

    private:
        mutable std::string access_name_;
        std::string scalartype_;
        std::string name_;
        std::string size_;
    };

    /**
     * @brief The vec_expr_infos class
     */
    class vec_expr_infos{
    public:

        vec_expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

        void add_infos(vec_infos const & infos){
            vecs_.push_back(infos);
        }

        std::list<vec_infos> const & vectors() const{
            return vecs_;
        }

        std::string generate() const{
            std::string res(expression_string_base_);
            for(std::list<vec_infos>::const_iterator it = vecs_.begin() ; it!= vecs_.end() ; ++it){
                replace_all_string(res,it->name(),it->access_name());
            }
            return res;
        }

    private:
        std::list<vec_infos> vecs_;
        mutable std::string expression_string_base_;
    };


    /**
     * @brief The scal_infos class
     */
    class scal_infos{
    public:
        scalar_infos(std::string const & scalartype,
                          std::string const & name,
                     bool is_passed_by_value) : scalartype_(scalartype), name_(name), is_passed_by_value_(is_passed_by_value)
        {
            if(is_passed_by_value) access_name_ = name_;
        }
        std::string const & scalartype() const{ return scalartype_; }
        std::string const & name() const{ return name_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string access_name() const { return access_name_; }

    private:
        mutable std::string access_name_;
        std::string scalartype_;
        std::string name_;
        bool is_passed_by_value_;
    };

    /**
     * @brief The scal_expr_infos class
     */
    class scal_expr_infos{
    public:

        scal_expr_infos(std::string const & expression_string_base) : expression_string_base_(expression_string_base){ }

        void add_infos(scal_infos const & infos){
            scals_.push_back(infos);
        }

        std::list<scal_infos> const & scalars() const{
            return scals_;
        }

        std::string generate() const{
            std::string res(expression_string_base_);
            for(std::list<scal_infos>::const_iterator it = scals_.begin() ; it!= scals_.end() ; ++it){
                replace_all_string(res,it->name(),it->access_name());
            }
            return res;
        }

    private:
        std::list<scal_infos> scals_;
        mutable std::string expression_string_base_;
    };

    /**
     * @brief The matvec_prod_infos struct
     */
    struct matvec_prod_infos{
        matvec_prod_infos(std::string const & name
                          ,vec_infos const & assigned
                          ,std::string const & op_str
                          ,mat_expr_infos const & lhs
                          ,vec_expr_infos const & rhs
                          ,vec_expr_infos const & additional_expression
                          ,std::string const & expression_string_base): name_(name)
                                                                        ,assigned_(assigned)
                                                                        ,lhs_(lhs)
                                                                         ,op_str_(op_str)
                                                                        , rhs_(rhs)
                                                                        , additional_expression_(additional_expression)
                                                                        , expression_string_base_(expression_string_base){ }
        mat_expr_infos const & lhs() const { return lhs_; }
        vec_expr_infos const & rhs() const { return rhs_; }
        vec_expr_infos const & additional_expression() const { return additional_expression_; }
        vec_infos const & assigned() const { return assigned_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string generate(){
            std::string prod(expression_string_base_);
            replace_all_string(prod,name_,access_name_);
            return assigned_.access_name() + op_str_ + prod + "+" + additional_expression_.generate();
        }

    private:
        std::string name_;
        std::string op_str_;
        vec_infos assigned_;
        mat_expr_infos lhs_;
        vec_expr_infos rhs_;
        vec_expr_infos additional_expression_;
        mutable std::string expression_string_base_;
        mutable std::string access_name_;
    };


    /**
     * @brief The inprod_infos struct
     */
    struct inprod_infos{
        inprod_infos(std::string const & name
                          ,scal_infos const & assigned
                          ,std::string const & op_str
                          ,vec_expr_infos const & lhs
                          ,vec_expr_infos const & rhs
                          ,scal_expr_infos const & additional_expression
                          ,std::string const & expression_string_base): name_(name)
                                                                        ,assigned_(assigned)
                                                                        ,lhs_(lhs)
                                                                         ,op_str_(op_str)
                                                                        , rhs_(rhs)
                                                                        , additional_expression_(additional_expression)
                                                                        , expression_string_base_(expression_string_base){ }
        vec_expr_infos const & lhs() const { return lhs_; }
        vec_expr_infos const & rhs() const { return rhs_; }
        scal_expr_infos const & additional_expression() const { return additional_expression_; }
        scal_infos const & assigned() const { return assigned_; }
        void access_name(std::string const & new_name) const{ access_name_ = new_name; }
        std::string generate(){
            std::string inprod(expression_string_base_);
            replace_all_string(prod,name_,access_name_);
            return assigned_.access_name() + op_str_ + inprod + "+" + additional_expression_.generate();
        }

    private:
        std::string name_;
        std::string op_str_;
        scal_infos assigned_;
        vec_expr_infos lhs_;
        vec_expr_infos rhs_;
        scal_expr_infos additional_expression_;
        mutable std::string expression_string_base_;
        mutable std::string access_name_;
    };

    template<class U>
    static vec_infos wrap_vector(){
         return vec_infos(print_type<typename U::ScalarType,1>::value(),
                                  U::name(),
                                  U::internal_size2_name());
    }

    /**
     * @brief wrap_vec_expr class
     */
    template<class T>
    struct wrap_vec_expr{
        typedef typename tree_utils::extract_if<T,result_of::is_symbolic_vector>::Result VectorsList;

        template<class U>
        struct functor{
            static void execute(vec_expr_infos & res){
                res.add_infos(wrap_vector<U>());
            }
        };

        static vec_expr_infos create(){
            vec_expr_infos res(make_expression_code<T>::value(""));
            typelist_utils::ForEach<VectorsList,functor>::execute(res);
            return res;
        }
    };

    template<class U>
    static mat_infos wrap_matrix(){
       return mat_infos(print_type<typename U::ScalarType,1>::value(),
                                      U::name(),
                                      U::size1_name(),
                                      U::size2_name(),
                                      result_of::is_row_major<U>::value,
                                      false);
    }

    template<class T>
    struct wrap_mat_expr{
        typedef typename tree_utils::extract_if<T,result_of::is_symbolic_matrix>::Result MatricesList;

        template<class U>
        struct functor{
            static void execute(mat_expr_infos & res){
                res.add_infos(wrap_matrix<U>());
            }
        };

        static mat_expr_infos create(){
            mat_expr_infos res(make_expression_code<T>::value(""));
            typelist_utils::ForEach<MatricesList,functor>::execute(res);
            return res;
        }
    };


    template<class T>
    struct wrap_matvec_prod{

        typedef typename tree_utils::remove_if<T,result_of::is_symbolic_vector,false >::Result Product;
        typedef typename tree_utils::remove_if<T,result_of::is_product_leaf,false >::Result Vectors;
        typedef typename Product::LHS ProdLHS;
        typedef typename Product::RHS ProdRHS;
        typedef typename result_of::expression_type<T>::Result    IntermediateType;
        static const unsigned int Alignment = IntermediateType::Alignment;

        static matvec_prod_infos create(){
            VIENNACL_STATIC_ASSERT(Alignment==1,AlignmentNotSupported);
            return matvec_prod_infos(Product::name()
                                     ,wrap_vector<typename T::LHS>()
                                     ,T::OP::expression_string()
                                     ,wrap_mat_expr<ProdLHS>::create()
                                     ,wrap_vec_expr<ProdRHS>::create()
                                     ,wrap_vec_expr<Vectors>::create()
                                     ,make_expression_code<Product>::value(""));
       }
    };


    std::string generate_matvec_prod_backend(std::vector<matvec_prod_infos> const & infos){
        std::string res;
        mat_infos const & mat = infos.front().lhs().matrices().front();
        vec_infos const & vec = infos.front().rhs().vectors().front();
        if (mat.is_rowmajor()){
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
//              generated_code+= U::private_value() + " += " + dot_product<LHS,RHS>::value("gid","gid") + ";\n";
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

