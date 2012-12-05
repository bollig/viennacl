#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_HPP

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


/** @file viennacl/generator/symbolic_types.hpp
 *  @brief Definition of the symbolic types.
 *
 *  Generator code contributed by Philippe Tillet
 */

//#include "viennacl/forwards.h"
//#include "viennacl/generator/result_of.hpp"
//#include "viennacl/generator/meta_tools/utils.hpp"

#include "viennacl/generator_fromscratch/utils.hpp"
#include "viennacl/forwards.h"

namespace viennacl
{
  namespace generator
  {

      template<class B,class T>
      class base_getter{
      public:
          static B & get(){
              static T res;
              return res;
          }
      };

      class infos_base{
      public:
          virtual std::string generate() const = 0;
          virtual std::string name() const = 0;
      };

      class binary_tree_infos_base{
      public:
          infos_base & lhs(){ return lhs_; }
          infos_base & rhs(){ return rhs_; }
      protected:
          binary_tree_infos_base(infos_base & lhs, infos_base & rhs) : lhs_(lhs), rhs_(rhs){        }
          infos_base & lhs_;
          infos_base & rhs_;
      };

      class unary_tree_infos_base{
      public:
          infos_base & sub(){ return sub_; }
      protected:
          unary_tree_infos_base(infos_base & sub) : sub_(sub){        }
          infos_base & sub_;
      };

      class leaf_infos_base : public infos_base{
      public:
          void access_name(std::string const & new_name) { access_name_ = new_name; }
          virtual ~leaf_infos_base(){ }
          bool is_modified(){ return is_modified_;}
          void is_modified(bool val){ is_modified_ = val; }
          virtual std::string generate() const { return access_name_; }
          virtual std::string name() const { return name_; }
      protected:
          leaf_infos_base(std::string const & scalartype,
                          std::string const & name): infos_base(), scalartype_(scalartype), name_(name), is_modified_(false){ }
          std::string scalartype_;
          std::string access_name_;
          std::string name_;
          bool is_modified_;
      };

      class kernel_argument : public leaf_infos_base{
      public:
          kernel_argument(int id, std::string const & scalartype, std::string const & name) : leaf_infos_base(scalartype,name), id_(id){ }
          int id() const{ return id_; }
          virtual std::string kernel_arguments() const = 0;
      private:
          int id_;
      };



      class arithmetic_tree_infos_base :  public infos_base,public binary_tree_infos_base{
      public:
          infos_base & op() { return op_; }
          std::string generate() const { return "(" + lhs_.generate() + op_.generate() + rhs_.generate() + ")"; }
          std::string name() const { return lhs_.name() + op_.name() + rhs_.name(); }
          arithmetic_tree_infos_base(infos_base & lhs, infos_base& op, infos_base & rhs) : binary_tree_infos_base(lhs,rhs), op_(op){        }
      private:
          infos_base & op_;
      };

      template<class LHS_, class OP_, class RHS_>
      class compound_node : public arithmetic_tree_infos_base, public base_getter<infos_base, compound_node<LHS_,OP_,RHS_> >
      {
        public:
          typedef LHS_  LHS;
          typedef RHS_  RHS;
          typedef OP_   OP;
          compound_node() : arithmetic_tree_infos_base(LHS::get(), OP::get(), RHS::get()){ }
      };

      template<class SUB_>
      class unary_minus : public unary_tree_infos_base, public virtual infos_base, public base_getter<infos_base, unary_minus<SUB_> >{
      public:
          unary_minus() : unary_tree_infos_base(SUB_::get()){ }
          std::string generate() const{
              return "(-" + sub_.generate() +")";
          }
      };

      class op_infos_base : public infos_base{
      public:
          std::string generate() const{ return expr_; }
          std::string name() const { return name_; }
      protected:
          op_infos_base(std::string const & expr, std::string const & name) : expr_(expr), name_(name){ }
      private:
          std::string expr_;
          std::string name_;
      };



      class assign_type : public op_infos_base, public base_getter<infos_base,assign_type>{
      public:
        assign_type() : op_infos_base(" = ", "eq"){ }
      };

      class add_type : public op_infos_base, public base_getter<infos_base,add_type>{
      public:
        add_type() : op_infos_base(" + ", "p"){ }
      };

      class inplace_add_type : public op_infos_base, public base_getter<infos_base,inplace_add_type>{
      public:
        inplace_add_type() : op_infos_base(" += ", "p_eq"){ }
      };

      class sub_type : public op_infos_base, public base_getter<infos_base,sub_type>{
      public:
        sub_type() : op_infos_base(" - ", "m"){ }
      };

      class inplace_sub_type : public op_infos_base, public base_getter<infos_base,inplace_sub_type>{
      public:
        inplace_sub_type() : op_infos_base(" -= ", "m_eq"){ }
      };

      class scal_mul_type : public op_infos_base, public base_getter<infos_base, scal_mul_type>{
      public:
        scal_mul_type() : op_infos_base(" * ", "mu"){ }
      };

      class inplace_scal_mul_type : public op_infos_base, public base_getter<infos_base, inplace_scal_mul_type>{
          inplace_scal_mul_type() : op_infos_base(" *= ", "mu_eq"){ }
      };


      class scal_div_type : public op_infos_base, public base_getter<infos_base,scal_div_type>{
        scal_div_type() : op_infos_base(" / ", "div"){ }
      };

      class inplace_scal_div_type :  public op_infos_base, public base_getter<infos_base,inplace_scal_div_type>{
          inplace_scal_div_type() : op_infos_base(" /= ", "div_eq"){ }
      };


      class elementwise_prod_type :  public op_infos_base, public base_getter<infos_base,elementwise_prod_type>{
          elementwise_prod_type() : op_infos_base(" * ", "ewp"){ }
      };

      class elementwise_div_type :  public op_infos_base, public base_getter<infos_base,elementwise_div_type>{
          elementwise_div_type() : op_infos_base(" / ", "ewd"){ }
      };



      class scal_infos_base : public kernel_argument{
      protected:
          scal_infos_base(std::string const & scalartype, std::string const & name, int id) : kernel_argument(id,scalartype,name){ }
      };

      class cpu_scal_infos_base : public scal_infos_base{
      protected:
          cpu_scal_infos_base(std::string const & scalartype, std::string const & name, int id) : scal_infos_base(scalartype,name,id){ }
      public:
          std::string kernel_arguments() const{
              return scalartype_ + " " + name_;
          }
      };

      class gpu_scal_infos_base : public scal_infos_base{
      protected:
          gpu_scal_infos_base(std::string const & scalartype, std::string const & name, int id) : scal_infos_base(scalartype,name,id){ }
      public:
          std::string kernel_arguments() const{
              return "__global " + scalartype_ + " " + name_;
          }
      };

      class inprod_infos_base : public scal_infos_base, public binary_tree_infos_base{
      public:
          enum step_t{reduce, compute};
          step_t step(){ return step_; }
      protected:
          inprod_infos_base(std::string const & scalartype, std::string const & name
                                 , infos_base & lhs, infos_base & rhs
                            ,step_t step) : scal_infos_base(scalartype,name,-1), binary_tree_infos_base(lhs,rhs), step_(step){        }
      private:
          step_t step_;
      };



//      template<class T>
//      class inprod_infos<T, typename viennacl::enable_if<result_of::is_inner_product_impl<T>::value>::type> : public inprod_infos_base{
//      typedef typename T::Type U;
//      public:
//          inprod_infos() : inprod_infos_base(print_type<typename U::ScalarType,1>::value(), T::name()
//                                                     ,vec_expr_infos<typename U::LHS>::get()
//                                                     , vec_expr_infos<typename U::RHS>::get()
//                                                     ,inprod_infos_base::compute){        }
//          static inprod_infos_base & get(){
//              static inprod_infos<T> res;
//              return res;
//          }
//      };

//      template<class T>
//      class inprod_infos<T, typename viennacl::enable_if<result_of::is_inner_product_leaf<T>::value>::type> : public inprod_infos_base{
//      public:
//          inprod_infos() : inprod_infos_base(print_type<typename T::ScalarType,1>::value("__local"), T::name()
//                                                     ,vec_expr_infos<typename T::LHS>::get()
//                                                     , vec_expr_infos<typename T::RHS>::get()
//                                                     ,inprod_infos_base::reduce){        }
//          static inprod_infos_base & get(){
//              static inprod_infos<T> res;
//              return res;
//          }
//      };

    /**
    * @brief Symbolic scalar type. Will be passed by value.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
    */
    template <unsigned int ID, typename SCALARTYPE>
    class cpu_symbolic_scalar : public cpu_scal_infos_base, public base_getter<infos_base, cpu_symbolic_scalar<ID,SCALARTYPE> >
    {
        cpu_symbolic_scalar() : cpu_scal_infos_base(print_type<SCALARTYPE,1>::value(),"c_s" + to_string(ID),ID){ }
    };

    /**
    * @brief Symbolic scalar type. Will be passed by pointer.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The SCALARTYPE of the scalar in the generated code
    */
    template <unsigned int ID, typename SCALARTYPE>
    class gpu_symbolic_scalar : public gpu_scal_infos_base, public base_getter<infos_base, gpu_symbolic_scalar<ID,SCALARTYPE> >
    {
      private:
        typedef gpu_symbolic_scalar<ID,SCALARTYPE> self_type;
        gpu_symbolic_scalar() : gpu_scal_infos_base(print_type<SCALARTYPE*,1>::value(), "g_s" + to_string(ID),ID ){ }

      public:
        typedef viennacl::scalar<SCALARTYPE> runtime_type;

        leaf_infos_base& get(){
            static self_type res;
            return res;
        }

        template<typename RHS_TYPE>
        compound_node<self_type, assign_type, RHS_TYPE >
        operator= ( RHS_TYPE const & ) const
        {
          return compound_node<self_type,assign_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        compound_node<self_type, inplace_scal_mul_type, RHS_TYPE >
        operator*= ( RHS_TYPE const & ) const
        {
          return compound_node<self_type,inplace_scal_mul_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        compound_node<self_type, inplace_scal_div_type, RHS_TYPE >
        operator/= ( RHS_TYPE const & ) const
        {
          return compound_node<self_type,inplace_scal_div_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        compound_node<self_type, inplace_add_type, RHS_TYPE >
        operator+= ( RHS_TYPE const & ) const
        {
          return compound_node<self_type,inplace_add_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        compound_node<self_type, inplace_sub_type, RHS_TYPE >
        operator-= ( RHS_TYPE const & ) const
        {
          return compound_node<self_type,inplace_sub_type,RHS_TYPE >();
        }

    };


      class vec_infos_base : public kernel_argument{
      public:
          std::string const & size() const{ return size_; }
          std::string const & internal_size() const{ return internal_size_;}
          std::string const & start() const{ return start_;}
          std::string const & inc() const{ return inc_;}
          std::string kernel_arguments() const
          {
            return " __global " + scalartype_ + " " + name_
                + ", unsigned int " + size_
                + ", unsigned int " + internal_size_+ "\n" ;
          }
          virtual ~vec_infos_base(){ }
      protected:
          vec_infos_base(std::string const & scalartype, std::string const & name, unsigned int id) :
                                                          kernel_argument(id,scalartype,name)
                                                          ,size_(name_+"_size"){ }
      private:
          std::string size_;
          std::string internal_size_;
          std::string start_;
          std::string inc_;
      };

      /**
      * @brief Symbolic vector type
      *
      * @tparam ID The argument ID of the vector in the generated code
      * @tparam SCALARTYPE The Scalartype of the vector in the generated code
      * @tparam ALIGNMENT The Alignment of the vector in the generated code
      */

      //TODO: Add start and inc...
      template <unsigned int ID, typename SCALARTYPE, unsigned int ALIGNMENT>
      class symbolic_vector : public vec_infos_base, public base_getter<infos_base, symbolic_vector<ID,SCALARTYPE,ALIGNMENT> >{
        private:
          typedef symbolic_vector<ID,SCALARTYPE,ALIGNMENT> self_type;

        public:

          typedef viennacl::vector<SCALARTYPE,ALIGNMENT> runtime_type;

          symbolic_vector() : vec_infos_base(print_type<SCALARTYPE,1>::value()
                                             , "v_a" + to_string(ALIGNMENT) + "_" + to_string(ID)
                                             , ID) { }


          template<typename RHS_TYPE>
          compound_node<self_type, assign_type, RHS_TYPE >
          operator= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,assign_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_scal_mul_type, RHS_TYPE >
          operator*= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_scal_mul_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_scal_div_type, RHS_TYPE >
          operator/= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_scal_div_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_add_type, RHS_TYPE >
          operator+= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_add_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_sub_type, RHS_TYPE >
          operator-= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_sub_type,RHS_TYPE >();
          }
      };




      class mat_infos_base : public kernel_argument{
      public:
          std::string const & size1() const{ return size1_; }
          std::string const & size2() const{ return size2_; }
          std::string const & internal_size1() const{ return internal_size1_; }
          std::string const & internal_size2() const{ return internal_size2_; }
          std::string const & row_inc() const{ return row_inc_; }
          std::string const & col_inc() const{ return col_inc_;}
          std::string const & row_start() const{ return row_start_;}
          std::string const & col_start() const{ return col_start_;}
          std::string kernel_arguments() const{
              return " __global " + scalartype_  + " " + name()
                    + ", unsigned int " + row_start_
                    + ", unsigned int " + col_start_
                    + ", unsigned int " + row_inc_
                    + ", unsigned int " + col_inc_
                    + ", unsigned int " + size1_
                    + ", unsigned int " + size2_
                    + ", unsigned int " + internal_size1_
                    + ", unsigned int " + internal_size2_
                    + "\n";
          }
          bool const is_rowmajor() const { return is_rowmajor_; }
          bool const is_transposed() const { return is_transposed_; }
          virtual ~mat_infos_base() { }
      protected:
          mat_infos_base(std::string const & scalartype
                         ,std::string const & name
                         ,bool is_rowmajor
                         ,bool is_transposed
                         ,unsigned int id) : kernel_argument(id, scalartype, name)
                                                ,size1_(name_ + "_size1")
                                                ,size2_(name_ + "size2")
                                                ,internal_size1_(size1_ + "_internal")
                                                ,internal_size2_(size2_ + "_internal")
                                                ,row_inc_(name_ + "_row_inc")
                                                ,col_inc_(name_+"_col_inc")
                                                ,row_start_(name_+"_row_start")
                                                ,col_start_(name_+"_col_start")
                                                ,is_rowmajor_(is_rowmajor)
                                                ,is_transposed_(is_transposed){ }
      private:
          std::string size1_;
          std::string size2_;
          std::string internal_size1_;
          std::string internal_size2_;
          std::string row_inc_;
          std::string col_inc_;
          std::string row_start_;
          std::string col_start_;
          bool is_rowmajor_;
          bool is_transposed_;
      };


      /**
      * @brief Symbolic matrix type
      *
      * @tparam ID The argument ID of the matrix in the generated code
      * @tparam SCALARTYPE The Scalartype of the matrix in the generated code
      * @tparam F The Layout of the matrix in the generated code
      * @tparam ALIGNMENT The Alignment of the matrix in the generated code
      */
      template<unsigned int ID, typename SCALARTYPE, class F, unsigned int ALIGNMENT>
      class symbolic_matrix : public mat_infos_base, public base_getter<infos_base, symbolic_matrix<ID, SCALARTYPE, F, ALIGNMENT> >
      {

          typedef symbolic_matrix<ID, SCALARTYPE, F, ALIGNMENT> self_type;

          static std::string name()
          {
            return "m_a_" + viennacl::generator::to_string(F()) + "_"
                          + viennacl::generator::to_string(ALIGNMENT) + "_"
                          + viennacl::generator::to_string<long>(ID);
          }


        public:
          symbolic_matrix() : mat_infos_base(print_type<SCALARTYPE,1>::value()
                                        ,name()
                                        ,true
                                        ,false){ }

          typedef viennacl::matrix<SCALARTYPE,F,ALIGNMENT> runtime_type;

          template<typename RHS_TYPE>
          compound_node<self_type, assign_type, RHS_TYPE >
          operator= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,assign_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_scal_mul_type, RHS_TYPE >
          operator*= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_scal_mul_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_scal_div_type, RHS_TYPE >
          operator/= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_scal_div_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_add_type, RHS_TYPE >
          operator+= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_add_type,RHS_TYPE >();
          }

          template<typename RHS_TYPE>
          compound_node<self_type, inplace_sub_type, RHS_TYPE >
          operator-= ( RHS_TYPE const & ) const
          {
            return compound_node<self_type,inplace_sub_type,RHS_TYPE >();
          }
      };


      /** @brief Unary minus operator */
      template<class T>
      unary_minus<T> operator -(T const &)
      {
        return unary_minus<T>();
      }

      /** @brief Scalar multiplication operator */
      template<class LHS_TYPE, class RHS_TYPE>
      compound_node<LHS_TYPE,scal_mul_type,RHS_TYPE>
      operator* ( LHS_TYPE const &, RHS_TYPE const & )
      {
        return compound_node<LHS_TYPE, scal_mul_type,RHS_TYPE> ();
      }



      /** @brief Scalar division operator */
      template<class LHS_TYPE, class RHS_TYPE>
      compound_node<LHS_TYPE,scal_div_type,RHS_TYPE>
      operator/ ( LHS_TYPE const &, RHS_TYPE const & )
      {
        return compound_node<LHS_TYPE,scal_div_type,RHS_TYPE> ();
      }

      /** @brief Addition operator on 2 elements of the same type */
      template<class LHS_TYPE, class RHS_TYPE>
      compound_node<LHS_TYPE, add_type, RHS_TYPE>
      operator+ ( LHS_TYPE const &, RHS_TYPE const & )
      {
        return compound_node<LHS_TYPE, add_type, RHS_TYPE>();
      }

      /** @brief Substraction operator on 2 elements of the same type */
      template<class LHS_TYPE, class RHS_TYPE>
      compound_node<LHS_TYPE, sub_type, RHS_TYPE>
      operator- ( LHS_TYPE const &, RHS_TYPE const & )
      {
        return compound_node<LHS_TYPE, sub_type, RHS_TYPE>();
      }


      template<class T>
      static bool is_type(infos_base* p){
          return dynamic_cast<T *>(p);
      }


//      template <class T>
//      struct inner_prod_impl_t
//      {
//        typedef T Type;

//        static std::string name()
//        {
//          return T::name();
//        }

//        static std::string private_value()
//        {
//          return "private_"+name();
//        }

//        static std::string declarations()
//        {
//          return print_type<typename T::ScalarType,1>::value() + " " + private_value() + "=0;\n" ;
//        }

//        static std::string kernel_arguments() const
//        {
//          return T::kernel_arguments();
//        }

//        enum { id = T::id };
//      };

//      template<class LHS_, class RHS_>
//      class compound_node<LHS_,inner_prod_type,RHS_>
//      {
//        public:
//          /**
//          * @brief Specialization for the inner product
//          */
//          typedef LHS_ LHS;
//          typedef RHS_ RHS;
//          typedef inner_prod_type OP;
//          typedef typename result_of::expression_type<LHS>::Result IntermediateType;  //Note: Visual Studio does not allow to combine this line with the next one directly.
//          typedef typename IntermediateType::ScalarType ScalarType;
//          static const unsigned int Alignment = IntermediateType::Alignment;

//          enum { id = -2 };

//          static std::string kernel_arguments() const
//          {
//              return  "__global " + print_type<ScalarType*,1>::value() + " " + name() + '\n';
//          }

//          static std::string name()
//          {
//              return  LHS::name() + "_inprod_" + RHS::name();
//          }

//          static std::string local_value()
//          {
//              return "local_"+name();
//          }

//          static std::string declarations()
//          {
//              return "__local " + print_type<ScalarType,1>::value() + " " + local_value() + ";\n";
//          }

//          static std::string scalar_name()
//          {
//            return name() +"_s";
//          }

//      };

//      /**
//      * @brief Specialization for the matrix-vector product.
//      */
//      template<class LHS_, class RHS_>
//      class compound_node<LHS_,prod_type,RHS_>
//      {
//        private:
//          typedef compound_node<LHS_,prod_type,RHS_> self_type;

//        public:
//          typedef LHS_ LHS;
//          typedef RHS_ RHS;

//          typedef prod_type OP;
//          enum { id = LHS::id };

//          typedef typename result_of::expression_type<LHS>::Result IntermediateType;    //Note: Visual Studio does not allow to combine this line with the next one directly.
//          typedef typename IntermediateType::ScalarType ScalarType;
//          static const unsigned int Alignment = result_of::expression_type<LHS>::Result::Alignment;

//          static std::string name()
//          {
//            return LHS::name() + "_prod_" + RHS::name();
//          }
//      };


  } // namespace generator
} // namespace viennacl


#endif
