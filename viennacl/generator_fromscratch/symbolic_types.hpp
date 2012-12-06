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
          virtual ~infos_base(){ }
      };

      class binary_tree_infos_base{
      public:
          infos_base & lhs(){ return *lhs_; }
          infos_base & rhs(){ return *rhs_; }
      protected:
          binary_tree_infos_base(infos_base * lhs, infos_base * rhs) : lhs_(lhs), rhs_(rhs){        }
          viennacl::tools::shared_ptr<infos_base> lhs_;
          viennacl::tools::shared_ptr<infos_base> rhs_;
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
          bool is_modified(){ return is_assigned_;}
          virtual std::string generate() const { return access_name_; }
          virtual std::string name() const { return name_; }
          std::string const & scalartype() const { return scalartype_; }
          bool operator<(leaf_infos_base const & other) const { return (&access_name_ < &other.access_name_); }

      protected:
          leaf_infos_base(std::string & access_name,
                          std::string const & scalartype,
                          std::string const & name): infos_base(), access_name_(access_name), scalartype_(scalartype),name_(name), is_assigned_(false){ }
          std::string & access_name_;
          std::string scalartype_;
          std::string name_;
          bool is_assigned_;
      };

      class kernel_argument : public leaf_infos_base{
      public:
          kernel_argument(std::string & access_name, std::string const & scalartype, std::string const & name,int id) : leaf_infos_base(access_name,scalartype,name), id_(id){ }
          int id() const{ return id_; }
          virtual std::string kernel_arguments() const = 0;
      private:
          int id_;
      };



      class arithmetic_tree_infos_base :  public infos_base,public binary_tree_infos_base{
      public:
          infos_base & op() { return op_; }
          std::string generate() const { return "(" + (*lhs_).generate() + op_.generate() + (*rhs_).generate() + ")"; }
          std::string name() const { return (*lhs_).name() + op_.name() + (*rhs_).name(); }
          arithmetic_tree_infos_base(infos_base * lhs, infos_base& op, infos_base * rhs) : binary_tree_infos_base(lhs,rhs), op_(op){        }
      private:
          infos_base & op_;
      };

      class vector_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          vector_expression_infos_base(infos_base * lhs, infos_base& op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };

      class scalar_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          scalar_expression_infos_base(infos_base * lhs, infos_base& op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };

      class matrix_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          matrix_expression_infos_base(infos_base * lhs, infos_base& op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };


      template<class LHS, class OP, class RHS>
      class vector_expression : public vector_expression_infos_base{
      public:
          vector_expression(LHS const & lhs, RHS const & rhs) :vector_expression_infos_base(new LHS(lhs),OP::get(),new RHS(rhs)){ }
      };

      template<class LHS, class OP, class RHS>
      class scalar_expression : public scalar_expression_infos_base{
      public:
          scalar_expression(LHS const & lhs, RHS const & rhs) :scalar_expression_infos_base(static_cast<infos_base &>(lhs_value_),OP::get(),static_cast<infos_base &>(rhs_value_))
            , lhs_value_(lhs), rhs_value_(rhs){ }
      private:
         LHS lhs_value_;
         RHS rhs_value_;
      };

      template<class LHS, class OP, class RHS>
      class matrix_expression : public matrix_expression_infos_base{
      public:
          matrix_expression(LHS const & lhs, RHS const & rhs) :matrix_expression_infos_base(static_cast<infos_base &>(lhs_value_),OP::get(),static_cast<infos_base &>(rhs_value_))
            , lhs_value_(lhs), rhs_value_(rhs){ }
      private:
         LHS lhs_value_;
         RHS rhs_value_;
      };



      template<class SUB_>
      class unary_minus : public unary_tree_infos_base, public infos_base{
      public:
          unary_minus(SUB_ sub) : unary_tree_infos_base(static_cast<infos_base &>(sub)){ }
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
          scal_infos_base(std::string & access_name,std::string const & scalartype, std::string const & name, int id) : kernel_argument(access_name,scalartype,name,id){ }
      };

      class cpu_scal_infos_base : public scal_infos_base{
      protected:
          cpu_scal_infos_base(std::string & access_name,std::string const & scalartype, std::string const & name, int id) : scal_infos_base(access_name,scalartype,name,id){ }
      public:
          std::string kernel_arguments() const{
              return scalartype_ + " " + name_;
          }
      };

      class gpu_scal_infos_base : public scal_infos_base{
      protected:
          gpu_scal_infos_base(std::string & access_name,std::string const & scalartype, std::string const & name, int id) : scal_infos_base(access_name,scalartype,name,id){ }
      public:
          std::string kernel_arguments() const{
              return "__global " + scalartype_ + " " + name_;
          }
      };

      class inprod_infos_base : public binary_tree_infos_base, public scal_infos_base{
      public:
          enum step_t{compute,reduce};
          step_t step(){ return step_; }
          void step(step_t s){ step_ = s; }

      protected:
          inprod_infos_base(infos_base * lhs, infos_base * rhs
                            ,std::string const & scalartype
                            ,std::string & access_name
                            ,step_t & step): binary_tree_infos_base(lhs,rhs)
                                             ,scal_infos_base(access_name
                                                              ,scalartype
                                                              ,lhs_->name() + "_inprod_" + rhs_->name()
                                                              ,-1)
                                            ,step_(step){ }
      private:
          step_t & step_;
      };

      template<class LHS, class RHS>
      class inprod_infos : public inprod_infos_base{
      public:

          inprod_infos(LHS const & lhs, RHS const & rhs) :
              inprod_infos_base(new LHS(lhs), new RHS(rhs),
                                print_type<typename LHS::ScalarType,1>::value(),
                                access_name_, step_){ }

          std::string kernel_arguments() const{
              return "__global " + scalartype_ + " " + name_;
          }
      private:
          static step_t step_;
          static std::string access_name_;
          LHS lhs_value_;
          RHS rhs_value_;
      };


      template <class LHS, class RHS>
      std::string inprod_infos<LHS,RHS>::access_name_;

      template <class LHS, class RHS>
      typename inprod_infos<LHS,RHS>::step_t inprod_infos<LHS,RHS>::step_ = inprod_infos<LHS,RHS>::compute;

    /**
    * @brief Symbolic scalar type. Will be passed by value.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
    */
    template <unsigned int ID, typename SCALARTYPE>
    class cpu_symbolic_scalar : public cpu_scal_infos_base
    {
    private:
        static std::string access_name_;
    public:
        cpu_symbolic_scalar() : cpu_scal_infos_base(access_name_,print_type<SCALARTYPE,1>::value(),"c_s" + to_string(ID),ID){ }
    };

    template <unsigned int ID, typename SCALARTYPE>
    std::string cpu_symbolic_scalar<ID,SCALARTYPE>::access_name_;

    /**
    * @brief Symbolic scalar type. Will be passed by pointer.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The SCALARTYPE of the scalar in the generated code
    */
    template <unsigned int ID, typename SCALARTYPE>
    class gpu_symbolic_scalar : public gpu_scal_infos_base
    {
      private:
        typedef gpu_symbolic_scalar<ID,SCALARTYPE> self_type;
        static std::string access_name_;

      public:
        typedef viennacl::scalar<SCALARTYPE> runtime_type;
        gpu_symbolic_scalar() : gpu_scal_infos_base(access_name_,print_type<SCALARTYPE*,1>::value(), "g_s" + to_string(ID),ID ){ }

        leaf_infos_base& get(){
            static self_type res;
            return res;
        }

        template<typename RHS_TYPE>
        scalar_expression<self_type, assign_type, RHS_TYPE >
        operator= ( RHS_TYPE const & rhs )
        {
          is_assigned_=true;
          return scalar_expression<self_type,assign_type,RHS_TYPE >(*this,rhs);
        }

        template<typename RHS_TYPE>
        scalar_expression<self_type, inplace_scal_mul_type, RHS_TYPE >
        operator*= ( RHS_TYPE const & rhs )
        {
          is_assigned_=true;
          return scalar_expression<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
        }

        template<typename RHS_TYPE>
        scalar_expression<self_type, inplace_scal_div_type, RHS_TYPE >
        operator/= ( RHS_TYPE const & rhs )
        {
          is_assigned_=true;
          return scalar_expression<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
        }

        template<typename RHS_TYPE>
        scalar_expression<self_type, inplace_add_type, RHS_TYPE >
        operator+= ( RHS_TYPE const & rhs )
        {
          is_assigned_=true;
          return scalar_expression<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
        }

        template<typename RHS_TYPE>
        scalar_expression<self_type, inplace_sub_type, RHS_TYPE >
        operator-= ( RHS_TYPE const & rhs )
        {
          is_assigned_=true;
          return scalar_expression<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
        }

    };

    template <unsigned int ID, typename SCALARTYPE>
    std::string gpu_symbolic_scalar<ID,SCALARTYPE>::access_name_;


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
          vec_infos_base(std::string & access_name, std::string const & scalartype, std::string const & name, unsigned int id) :
                                                          kernel_argument(access_name,scalartype,name,id)
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
      template <unsigned int ID, typename SCALARTYPE, unsigned int ALIGNMENT=1>
      class symbolic_vector : public vec_infos_base{
        private:
          typedef symbolic_vector<ID,SCALARTYPE,ALIGNMENT> self_type;
          static std::string access_name_;
        public:

          typedef viennacl::vector<SCALARTYPE,ALIGNMENT> runtime_type;

          symbolic_vector() : vec_infos_base(access_name_,
                                             print_type<SCALARTYPE,1>::value()
                                             , "v_a" + to_string(ALIGNMENT) + "_" + to_string(ID)
                                             , ID) { }


          template<typename RHS_TYPE>
          vector_expression<self_type, assign_type, RHS_TYPE >
          operator= ( RHS_TYPE const & rhs)
          {
            is_assigned_=true;
            return vector_expression<self_type,assign_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          vector_expression<self_type, inplace_scal_mul_type, RHS_TYPE >
          operator*= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return vector_expression<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          vector_expression<self_type, inplace_scal_div_type, RHS_TYPE >
          operator/= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return vector_expression<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          vector_expression<self_type, inplace_add_type, RHS_TYPE >
          operator+= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return vector_expression<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          vector_expression<self_type, inplace_sub_type, RHS_TYPE >
          operator-= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return vector_expression<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
          }
      };

      template <unsigned int ID, typename SCALARTYPE, unsigned int ALIGNMENT>
      std::string symbolic_vector<ID,SCALARTYPE,ALIGNMENT>::access_name_;


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
          mat_infos_base(std::string & access_name,
                         std::string const & scalartype
                         ,std::string const & name
                         ,bool is_rowmajor
                         ,bool is_transposed
                         ,unsigned int id) : kernel_argument(access_name, scalartype, name,id)
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
      template<unsigned int ID, typename SCALARTYPE, class F, unsigned int ALIGNMENT=1>
      class symbolic_matrix : public mat_infos_base
      {

          typedef symbolic_matrix<ID, SCALARTYPE, F, ALIGNMENT> self_type;

          static std::string name()
          {
            return "m_a_" + viennacl::generator::to_string(F()) + "_"
                          + viennacl::generator::to_string(ALIGNMENT) + "_"
                          + viennacl::generator::to_string<long>(ID);
          }


          static std::string access_name_;
        public:
          symbolic_matrix() : mat_infos_base(access_name_,
                                             print_type<SCALARTYPE,1>::value()
                                        ,name()
                                        ,true
                                        ,false){ }

          typedef viennacl::matrix<SCALARTYPE,F,ALIGNMENT> runtime_type;

          template<typename RHS_TYPE>
          matrix_expression<self_type, assign_type, RHS_TYPE >
          operator= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return matrix_expression<self_type,assign_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          matrix_expression<self_type, inplace_scal_mul_type, RHS_TYPE >
          operator*= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return matrix_expression<self_type,inplace_scal_mul_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          matrix_expression<self_type, inplace_scal_div_type, RHS_TYPE >
          operator/= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return matrix_expression<self_type,inplace_scal_div_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          matrix_expression<self_type, inplace_add_type, RHS_TYPE >
          operator+= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return matrix_expression<self_type,inplace_add_type,RHS_TYPE >(*this,rhs);
          }

          template<typename RHS_TYPE>
          matrix_expression<self_type, inplace_sub_type, RHS_TYPE >
          operator-= ( RHS_TYPE const & rhs )
          {
            is_assigned_=true;
            return matrix_expression<self_type,inplace_sub_type,RHS_TYPE >(*this,rhs);
          }
      };

      template <unsigned int ID, typename SCALARTYPE, class F, unsigned int ALIGNMENT>
      std::string symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT>::access_name_;

      template<class T>
      struct is_vector_expression{ enum { value = 0 }; };

      template <unsigned int ID, typename SCALARTYPE, unsigned int ALIGNMENT>
      struct is_vector_expression<symbolic_vector<ID,SCALARTYPE,ALIGNMENT> >{ enum { value = 1 }; };

      template<class LHS, class OP, class RHS>
      struct is_vector_expression<vector_expression<LHS,OP,RHS> >{ enum { value = 1 }; };


      template<class T>
      struct is_matrix_expression{ enum { value = 0 }; };

      template <unsigned int ID, typename SCALARTYPE, class F,unsigned int ALIGNMENT>
      struct is_matrix_expression<symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> >{ enum { value = 1 }; };

      template<class LHS, class OP, class RHS>
      struct is_matrix_expression<matrix_expression<LHS,OP,RHS> >{ enum { value = 1 }; };

      template<class T>
      struct is_scalar_expression{ enum { value = 0 }; };

      template <unsigned int ID, typename SCALARTYPE>
      struct is_scalar_expression<gpu_symbolic_scalar<ID,SCALARTYPE> >{ enum { value = 1 }; };

      template <unsigned int ID, typename SCALARTYPE>
      struct is_scalar_expression<cpu_symbolic_scalar<ID,SCALARTYPE> >{ enum { value = 1 }; };

      template<class LHS, class OP, class RHS>
      struct is_scalar_expression<scalar_expression<LHS,OP,RHS> >{ enum { value = 1 }; };

      template<class LHS, class RHS>
      struct is_scalar_expression<inprod_infos<LHS,RHS> >{ enum { value = 1 }; };




      /** @brief Unary minus operator */
      template<class T>
      unary_minus<T> operator -(T const &)
      {
        return unary_minus<T>();
      }

//      /** @brief Scalar multiplication operator */
//      template<class LHS_TYPE, class RHS_TYPE>
//      compound_node<LHS_TYPE,scal_mul_type,RHS_TYPE>
//      operator* ( LHS_TYPE const &, RHS_TYPE const & )
//      {
//        return compound_node<LHS_TYPE, scal_mul_type,RHS_TYPE> ();
//      }



//      /** @brief Scalar division operator */
//      template<class LHS_TYPE, class RHS_TYPE>
//      compound_node<LHS_TYPE,scal_div_type,RHS_TYPE>
//      operator/ ( LHS_TYPE const &, RHS_TYPE const & )
//      {
//        return compound_node<LHS_TYPE,scal_div_type,RHS_TYPE> ();
//      }

      /**
       * Operator* for Scalar*Vector
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_vector_expression<RHS_TYPE>::value)
                                 ||(is_scalar_expression<RHS_TYPE>::value && is_vector_expression<LHS_TYPE>::value)
                                  ,vector_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator* ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return vector_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }

      /**
       * Operator* for Scalar*Matrix
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_matrix_expression<RHS_TYPE>::value)
                                 ||(is_scalar_expression<RHS_TYPE>::value && is_matrix_expression<LHS_TYPE>::value)
                                  ,matrix_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator* ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return matrix_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }

      /**
       * Operator* for Scalar*Scalar
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_scalar_expression<RHS_TYPE>::value)
                                  ,scalar_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator* ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return scalar_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * operator/ for Scalar*Vector
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_vector_expression<RHS_TYPE>::value)
                                 ||(is_scalar_expression<RHS_TYPE>::value && is_vector_expression<LHS_TYPE>::value)
                                  ,vector_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator/ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return vector_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }

      /**
       * operator/ for Scalar*Matrix
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_matrix_expression<RHS_TYPE>::value)
                                 ||(is_scalar_expression<RHS_TYPE>::value && is_matrix_expression<LHS_TYPE>::value)
                                  ,matrix_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator/ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return matrix_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }

      /**
       * operator/ for Scalar*Scalar
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<(is_scalar_expression<LHS_TYPE>::value && is_scalar_expression<RHS_TYPE>::value)
                                  ,scalar_expression<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
      operator/ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return scalar_expression<LHS_TYPE,scal_mul_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * Operator+ for adding 2vectors
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_vector_expression<LHS_TYPE>::value && is_vector_expression<RHS_TYPE>::value
                                  ,vector_expression<LHS_TYPE,add_type,RHS_TYPE> >::type
      operator+ ( LHS_TYPE const & lhs , RHS_TYPE const & rhs){
          return vector_expression<LHS_TYPE,add_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * Operator+ for adding 2matrices
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_matrix_expression<LHS_TYPE>::value && is_matrix_expression<RHS_TYPE>::value
                                  ,matrix_expression<LHS_TYPE,add_type,RHS_TYPE> >::type
      operator+ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return matrix_expression<LHS_TYPE,add_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * Operator+ for adding 2scalars
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_scalar_expression<LHS_TYPE>::value && is_scalar_expression<RHS_TYPE>::value
                                  ,scalar_expression<LHS_TYPE,add_type,RHS_TYPE> >::type
      operator+ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return scalar_expression<LHS_TYPE,add_type,RHS_TYPE>(lhs,rhs);
      }




      /**
       * Operator- for adding 2vectors
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_vector_expression<LHS_TYPE>::value && is_vector_expression<RHS_TYPE>::value
                                  ,vector_expression<LHS_TYPE,sub_type,RHS_TYPE> >::type
      operator- ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return vector_expression<LHS_TYPE,sub_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * Operator- for adding 2matrices
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_matrix_expression<LHS_TYPE>::value && is_matrix_expression<RHS_TYPE>::value
                                  ,matrix_expression<LHS_TYPE,sub_type,RHS_TYPE> >::type
      operator- ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return matrix_expression<LHS_TYPE,sub_type,RHS_TYPE>(lhs,rhs);
      }


      /**
       * Operator- for adding 2scalars
       */
      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_scalar_expression<LHS_TYPE>::value && is_scalar_expression<RHS_TYPE>::value
                                  ,scalar_expression<LHS_TYPE,sub_type,RHS_TYPE> >::type
      operator- ( LHS_TYPE const &, RHS_TYPE const & ){
          return scalar_expression<LHS_TYPE,sub_type,RHS_TYPE>();
      }



      template<class LHS_TYPE, class RHS_TYPE>
      typename viennacl::enable_if<is_vector_expression<LHS_TYPE>::value && is_vector_expression<RHS_TYPE>::value
                                  ,inprod_infos<LHS_TYPE,RHS_TYPE> >::type
      inner_prod ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ){
          return inprod_infos<LHS_TYPE,RHS_TYPE>(lhs,rhs);
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
