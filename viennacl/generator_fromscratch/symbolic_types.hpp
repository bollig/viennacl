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



#include "viennacl/generator_fromscratch/dummy_types.hpp"
#include "viennacl/generator_fromscratch/symbolic_types_base.hpp"

namespace viennacl
{
  namespace generator
  {






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

      template<class LHS, class OP, class RHS>
      class vector_expression : public vector_expression_infos_base{
      public:
          vector_expression(LHS const & lhs, RHS const & rhs) :vector_expression_infos_base( new LHS(lhs),new OP(),new RHS(rhs)){ }
      };

      template<class LHS, class OP, class RHS>
      class scalar_expression : public scalar_expression_infos_base{
      public:
          scalar_expression(LHS const & lhs, RHS const & rhs) :scalar_expression_infos_base( new LHS(lhs),new OP(),new RHS(rhs)){ }
      };

      template<class LHS, class OP, class RHS>
      class matrix_expression : public matrix_expression_infos_base{
      public:
          matrix_expression(LHS const & lhs, RHS const & rhs) :matrix_expression_infos_base( new LHS(lhs),new OP(),new RHS(rhs)){ }
      };



      template<class SUB_>
      class unary_minus : public infos_base, public unary_tree_infos_base{
      public:
          unary_minus(SUB_ sub) : unary_tree_infos_base(new SUB_(sub)){ }

      };



      template<class LHS, class RHS>
      class inprod_infos : public inprod_infos_base{
      public:
          inprod_infos(LHS const & lhs, RHS const & rhs) :
              inprod_infos_base(
                                new LHS(lhs), new RHS(rhs),new step_t(inprod_infos_base::compute),
                                print_type<typename LHS::ScalarType>::value()){ }

          std::string kernel_arguments() const{
              return "__global " + scalartype_ + "*" + " " + name_;
          }
      private:
          viennacl::vector<typename LHS::ScalarType> temp_;
      };

    /**
    * @brief Symbolic scalar type. Will be passed by value.
    *
    * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
    */
    template <typename SCALARTYPE>
    class cpu_symbolic_scalar : public cpu_scal_infos_base
    {
    public:
        typedef SCALARTYPE ScalarType;
        cpu_symbolic_scalar(ScalarType val) : cpu_scal_infos_base(print_type<SCALARTYPE>::value()), val_(val){ }
    private:
        ScalarType val_;
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
        gpu_symbolic_scalar() : gpu_scal_infos_base(print_type<SCALARTYPE>::value()){ }

        leaf_infos_base& get(){
            static self_type res;
            return res;
        }

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
          symbolic_vector(std::map<viennacl::backend::mem_handle, std::string> & access_names_map,
                          runtime_type const & rt) : vec_infos_base( print_type<SCALARTYPE>::value()), rt_obj_(rt){
              access_name_ = &access_names_map.insert(std::make_pair(rt_obj_.handle(),std::string())).first->second;
              name_ = "arg_"+to_string(access_names_map.size());
          }
        private:
          runtime_type const & rt_obj_;
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
          symbolic_matrix() : mat_infos_base( print_type<SCALARTYPE>::value()
                                        ,true
                                        ,false){ }

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
      struct dummy2exptree_impl
      {
      private:
        typedef typename T::Lhs Lhs;
        typedef typename T::Rhs Rhs;
        typedef typename dummy2exptree_impl<Lhs>::result_type LhsResult;
        typedef typename dummy2exptree_impl<Rhs>::result_type RhsResult;
      public:
        typedef typename get_symbolic_type<T,LhsResult,RhsResult>::type result_type;

          static result_type execute(std::map<viennacl::backend::mem_handle,std::string> & access_names_map, T const & t){
              return result_type(dummy2exptree_impl<Lhs>::execute(access_names_map, t.lhs())
                                 ,dummy2exptree_impl<Rhs>::execute(access_names_map, t.rhs()));
          }
      };

      template<class ScalarType, unsigned int Alignment>
      struct dummy2exptree_impl<dummy_vector<ScalarType, Alignment> >{
          typedef symbolic_vector<ScalarType, Alignment> result_type;
          static result_type execute(std::map<viennacl::backend::mem_handle,std::string> & access_names_map, dummy_vector<ScalarType,Alignment> const & v){
              return result_type(access_names_map,v.vec());
          }
      };

      template<class T>
      typename dummy2exptree_impl<T>::result_type dummy2exptree(std::map<viennacl::backend::mem_handle,std::string> & access_names_map, T const & t){
          return dummy2exptree_impl<T>::execute(access_names_map,t);
      }




  } // namespace generator
} // namespace viennacl


#endif
