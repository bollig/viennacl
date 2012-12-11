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
#include "viennacl/generator_fromscratch/dummy_types.hpp"

namespace viennacl
{
  namespace generator
  {

      class local_memory{
      public:

          local_memory(std::string const & name, unsigned int size, std::string const & scalartype): name_(name), size_(size), scalartype_(scalartype){ }

          std::string declare() const{
              return "__local " + scalartype_ + " " + name_ + '[' + to_string(size_) + ']';
          }

          unsigned int size() const{ return size_; }

          std::string access(std::string const & index) const{
              return name_ + '[' + index + ']';
          }

      private:
          std::string name_;
          unsigned int size_;
          std::string const & scalartype_;
      };

      class infos_base{
      public:
          virtual std::string generate() const = 0;
          virtual std::string name() const = 0;
          virtual ~infos_base(){ }
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



      class assign_type : public op_infos_base{
      public:
        assign_type() : op_infos_base(" = ", "eq"){ }
      };

      class add_type : public op_infos_base{
      public:
        add_type() : op_infos_base(" + ", "p"){ }
      };

      class inplace_add_type : public op_infos_base{
      public:
        inplace_add_type() : op_infos_base(" += ", "p_eq"){ }
      };

      class sub_type : public op_infos_base{
      public:
        sub_type() : op_infos_base(" - ", "m"){ }
      };

      class inplace_sub_type : public op_infos_base{
      public:
        inplace_sub_type() : op_infos_base(" -= ", "m_eq"){ }
      };

      class scal_mul_type : public op_infos_base{
      public:
        scal_mul_type() : op_infos_base(" * ", "mu"){ }
      };

      class inplace_scal_mul_type : public op_infos_base{
          inplace_scal_mul_type() : op_infos_base(" *= ", "mu_eq"){ }
      };


      class scal_div_type : public op_infos_base{
        scal_div_type() : op_infos_base(" / ", "div"){ }
      };

      class inplace_scal_div_type :  public op_infos_base{
          inplace_scal_div_type() : op_infos_base(" /= ", "div_eq"){ }
      };


      class elementwise_prod_type :  public op_infos_base{
          elementwise_prod_type() : op_infos_base(" * ", "ewp"){ }
      };

      class elementwise_div_type :  public op_infos_base{
          elementwise_div_type() : op_infos_base(" / ", "ewd"){ }
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
          void access_name(std::string const & new_name) { *access_name_ = new_name; }
          virtual ~leaf_infos_base(){ }
          bool is_modified(){ return is_assigned_;}
          virtual std::string generate() const { return *access_name_; }
          virtual std::string name() const { return name_; }
          std::string const & scalartype() const { return scalartype_; }
          bool operator<(leaf_infos_base const & other) const { return (access_name_ < other.access_name_); }

      protected:
          leaf_infos_base(std::string const & scalartype): infos_base(), scalartype_(scalartype){ }
          std::string * access_name_;
          std::string scalartype_;
          std::string name_;
          bool is_assigned_;
      };

      class kernel_argument : public leaf_infos_base{
      public:
          kernel_argument(std::string const & scalartype, int id) : leaf_infos_base(scalartype), id_(id){
              access_name_ = &access_names_map_.insert(std::make_pair(id,std::string())).first->second;
              name_ = "arg"+to_string(id_);
          }
          int id() const{ return id_; }
          virtual std::string kernel_arguments() const = 0;
      private:
          int id_;
          static std::map<unsigned int, std::string> access_names_map_;
      };

      std::map<unsigned int, std::string> kernel_argument::access_names_map_;


      class arithmetic_tree_infos_base :  public infos_base,public binary_tree_infos_base{
      public:
          infos_base & op() { return *op_; }
          std::string generate() const { return "(" + lhs_->generate() + op_->generate() + rhs_->generate() + ")"; }
          std::string name() const { return lhs_->name() + op_->name() + rhs_->name(); }
          arithmetic_tree_infos_base(infos_base * lhs, infos_base* op, infos_base * rhs) : binary_tree_infos_base(lhs,rhs), op_(op){        }
      private:
          viennacl::tools::shared_ptr<infos_base> op_;
      };

      class vector_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          vector_expression_infos_base(infos_base * lhs, infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };

      class scalar_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          scalar_expression_infos_base(infos_base * lhs, infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };

      class matrix_expression_infos_base : public arithmetic_tree_infos_base{
      public:
          matrix_expression_infos_base(infos_base * lhs, infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base(lhs,op,rhs){ }
      };


      template<class LHS, class OP, class RHS>
      class vector_expression : public vector_expression_infos_base{
      public:
          vector_expression(LHS const & lhs, RHS const & rhs) :vector_expression_infos_base(new LHS(lhs),new OP(),new RHS(rhs)){ }
      };

      template<class LHS, class OP, class RHS>
      class scalar_expression : public scalar_expression_infos_base{
      public:
          scalar_expression(LHS const & lhs, RHS const & rhs) :scalar_expression_infos_base(new LHS(lhs),new OP(),new RHS(rhs)){ }
      };

      template<class LHS, class OP, class RHS>
      class matrix_expression : public matrix_expression_infos_base{
      public:
          matrix_expression(LHS const & lhs, RHS const & rhs) :matrix_expression_infos_base(new LHS(lhs),new OP(),new RHS(rhs)){ }
      };



      template<class SUB_>
      class unary_minus : public unary_tree_infos_base, public infos_base{
      public:
          unary_minus(SUB_ sub) : unary_tree_infos_base(static_cast<infos_base &>(sub)){ }
          std::string generate() const{
              return "(-" + sub_.generate() +")";
          }
      };


      class scal_infos_base : public kernel_argument{
      protected:
          scal_infos_base(std::string const & scalartype, int id) : kernel_argument(scalartype,id){ }
      };

      class cpu_scal_infos_base : public scal_infos_base{
      protected:
          cpu_scal_infos_base(std::string const & scalartype
                              , int id) : scal_infos_base(scalartype,id){ }
      public:
          std::string kernel_arguments() const{
              return scalartype_ + " " + name_;
          }
      };

      class gpu_scal_infos_base : public scal_infos_base{
      protected:
          gpu_scal_infos_base(std::string const & scalartype, int id) : scal_infos_base(scalartype,id){ }
      public:
          std::string kernel_arguments() const{
              return "__global " + scalartype_ + "*"  + " " + name_;
          }
      };

      class inprod_infos_base : public binary_tree_infos_base, public scal_infos_base{
      public:
          enum step_t{compute,reduce};
          step_t step(){ return step_; }
          void step(step_t s){ step_ = s; }
          local_memory make_local_memory(unsigned int size){
              return local_memory(name_,size,scalartype_);
          }

      protected:
          inprod_infos_base(infos_base * lhs, infos_base * rhs
                            ,std::string const & scalartype
                            ,step_t & step): binary_tree_infos_base(lhs,rhs)
                                             ,scal_infos_base(scalartype,-1)
                                            ,step_(step){ }


      private:
          step_t & step_;
      };

      template<class LHS, class RHS>
      class inprod_infos : public inprod_infos_base{
      public:

          inprod_infos(LHS const & lhs, RHS const & rhs) :
              inprod_infos_base(new LHS(lhs), new RHS(rhs),
                                print_type<typename LHS::ScalarType>::value(), step_){ }

          std::string kernel_arguments() const{
              return "__global " + scalartype_ + "*" + " " + name_;
          }
      private:
          static step_t step_;
      };

      template <class LHS, class RHS>
      typename inprod_infos<LHS,RHS>::step_t inprod_infos<LHS,RHS>::step_ = inprod_infos<LHS,RHS>::compute;

    /**
    * @brief Symbolic scalar type. Will be passed by value.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
    */
    template <typename SCALARTYPE>
    class cpu_symbolic_scalar : public cpu_scal_infos_base
    {
    public:
        typedef SCALARTYPE ScalarType;
        cpu_symbolic_scalar(unsigned int id) : cpu_scal_infos_base(print_type<SCALARTYPE>::value()
                                                    ,"c_s" + to_string(id)
                                                    ,id){ }
    };


    /**
    * @brief Symbolic scalar type. Will be passed by pointer.
    *
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The SCALARTYPE of the scalar in the generated code
    */
    template <typename SCALARTYPE>
    class gpu_symbolic_scalar : public gpu_scal_infos_base
    {
      private:
        typedef gpu_symbolic_scalar<SCALARTYPE> self_type;

      public:
        typedef viennacl::scalar<SCALARTYPE> runtime_type;
        typedef SCALARTYPE ScalarType;
        gpu_symbolic_scalar(unsigned int id) : gpu_scal_infos_base(print_type<SCALARTYPE>::value(),id ){ }

        leaf_infos_base& get(){
            static self_type res;
            return res;
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
            return " __global " + scalartype_ + "*"  + " " + name_
                + ", unsigned int " + size_
                + ", unsigned int " + internal_size_+ "\n" ;
          }
          virtual ~vec_infos_base(){ }
      protected:
          vec_infos_base(std::string const & scalartype, unsigned int id) :
                                                          kernel_argument(scalartype,id)
                                                        ,size_(name_+"_size")
                                                        ,internal_size_(""){ }
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
      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      class symbolic_vector : public vec_infos_base{
        private:
          typedef symbolic_vector<SCALARTYPE,ALIGNMENT> self_type;
        public:
          typedef viennacl::vector<SCALARTYPE,ALIGNMENT> runtime_type;
          typedef SCALARTYPE ScalarType;
          symbolic_vector(runtime_type const & rt) : vec_infos_base(print_type<SCALARTYPE>::value(), 0) { }
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
              return " __global " + scalartype_ + "*"  + " " + name_
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
                         ,bool is_rowmajor
                         ,bool is_transposed
                         ,unsigned int id) : kernel_argument(scalartype,id)
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
      template<typename SCALARTYPE, class F, unsigned int ALIGNMENT=1>
      class symbolic_matrix : public mat_infos_base
      {
          typedef symbolic_matrix<SCALARTYPE, F, ALIGNMENT> self_type;

        public:
          symbolic_matrix(unsigned int id) : mat_infos_base(print_type<SCALARTYPE>::value()
                                        ,true
                                        ,false
                                        ,id){ }

          typedef viennacl::matrix<SCALARTYPE,F,ALIGNMENT> runtime_type;
          typedef SCALARTYPE ScalarType;
      };

      template<class Model, class LHS, class RHS>
      struct get_symbolic_type;

      template<class OP, class LHS1, class RHS1, class LHS2, class RHS2>
      struct get_symbolic_type<vector_expression_wrapper<LHS1,OP,RHS1>,LHS2,RHS2 >{ typedef vector_expression<LHS2,OP,RHS2> type;};
      template<class OP, class LHS1, class RHS1, class LHS2, class RHS2>
      struct get_symbolic_type<scalar_expression_wrapper<LHS1,OP,RHS1>,LHS2,RHS2 >{ typedef scalar_expression<LHS2,OP,RHS2> type;};
      template<class OP, class LHS1, class RHS1, class LHS2, class RHS2>
      struct get_symbolic_type<matrix_expression_wrapper<LHS1,OP,RHS1>,LHS2,RHS2 >{ typedef matrix_expression<LHS2,OP,RHS2> type;};

      template<class T>
      struct dummy2exptree_impl;

      template<class T>
      struct dummy2exptree_impl
      {
      private:
        typedef typename T::Lhs Lhs;
        typedef typename T::Rhs Rhs;
        typedef typename dummy2exptree_impl<Lhs>::result_type LhsResult;
        typedef typename dummy2exptree_impl<Rhs>::result_type RhsResult;
      public:
        typedef typename get_symbolic_type<T,LhsResult,RhsResult>::type result_type;

          static result_type execute(T const & t){
              return result_type(dummy2exptree_impl<Lhs>::execute(t.lhs())
                                 ,dummy2exptree_impl<Rhs>::execute(t.rhs()));
          }
      };

      template<class ScalarType, unsigned int Alignment>
      struct dummy2exptree_impl<dummy_vector<ScalarType, Alignment> >{
          typedef symbolic_vector<ScalarType, Alignment> result_type;
          static result_type execute(dummy_vector<ScalarType,Alignment> const & v){
              return result_type(v.vec());
          }
      };

      template<class T>
      typename dummy2exptree_impl<T>::result_type dummy2exptree(T const & t){
          return dummy2exptree_impl<T>::execute(t);
      }




  } // namespace generator
} // namespace viennacl


#endif
