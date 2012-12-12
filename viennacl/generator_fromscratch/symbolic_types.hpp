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
#include "viennacl/generator_fromscratch/functions.hpp"

namespace viennacl
{
  namespace generator
  {


      struct shared_infos{
      public:
          shared_infos(std::string access_name, std::string name) : access_name_(name), name_(name){ }
          std::string & access_name(){ return access_name_; }
          std::string const & name() const{ return name_; }
      private:
          std::string access_name_;
          std::string name_;
      };


      class assign_type : public op_infos_base{
      public:
          assign_type() : op_infos_base(" = ", "eq",true){ }
      };

      class add_type : public op_infos_base{
      public:
        add_type() : op_infos_base(" + ", "p",false){ }
      };

      class inplace_add_type : public op_infos_base{
      public:
        inplace_add_type() : op_infos_base(" += ", "p_eq",true){ }
      };

      class sub_type : public op_infos_base{
      public:
        sub_type() : op_infos_base(" - ", "m",false){ }
      };

      class inplace_sub_type : public op_infos_base{
      public:
        inplace_sub_type() : op_infos_base(" -= ", "m_eq",true){ }
      };

      class scal_mul_type : public op_infos_base{
      public:
        scal_mul_type() : op_infos_base(" * ", "mu",false){ }
      };

      class inplace_scal_mul_type : public op_infos_base{
          inplace_scal_mul_type() : op_infos_base(" *= ", "mu_eq",true){ }
      };


      class scal_div_type : public op_infos_base{
        scal_div_type() : op_infos_base(" / ", "div", false){ }
      };

      class inplace_scal_div_type :  public op_infos_base{
          inplace_scal_div_type() : op_infos_base(" /= ", "div_eq", true){ }
      };


      class elementwise_prod_type :  public op_infos_base{
          elementwise_prod_type() : op_infos_base(" * ", "ewp",false){ }
      };

      class elementwise_div_type :  public op_infos_base{
          elementwise_div_type() : op_infos_base(" / ", "ewd", false){ }
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
          symbolic_vector(std::string & access_name
                          ,std::string const & name
                          ,runtime_type const & rt) : vec_infos_base(name,print_type<SCALARTYPE>::value()), rt_obj_(rt){
              access_name_ = &access_name;
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

          static result_type execute(std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map,
                                     T const & t){
              return result_type(dummy2exptree_impl<Lhs>::execute(access_names_map, t.lhs())
                                 ,dummy2exptree_impl<Rhs>::execute(access_names_map, t.rhs()));
          }
      };

      template<class ScalarType, unsigned int Alignment>
      struct dummy2exptree_impl<dummy_vector<ScalarType, Alignment> >{
          typedef symbolic_vector<ScalarType, Alignment> result_type;
          static result_type execute(std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map,
                                     dummy_vector<ScalarType,Alignment> const & v){
              std::map<viennacl::backend::mem_handle,shared_infos>::iterator it = access_names_map.insert(std::make_pair(v.vec().handle(),shared_infos(std::string(), std::string("arg"+to_string(access_names_map.size()))))).first;
              return result_type(it->second.access_name(), it->second.name(), v.vec());
          }
      };




      template<class T1, class T2, class T3, class T4, class T5>
      struct dummy2exptree_impl<function_wrapper_impl<T1,T2,T3,T4,T5> >{
      private:
          template<class U>
          static void handle_function_arg(symbolic_function & fun, U const * t, std::string name
                              , std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map)
          {
              fun.add_arg(name, dummy2exptree_impl<U>::execute(access_names_map,*t));
          }

          static void handle_function_arg(symbolic_function & fun, void const* t, std::string name
                              , std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map)
          { }

      public:
          typedef symbolic_function result_type;
          static result_type execute(std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map
                                     ,function_wrapper_impl<T1,T2,T3,T4,T5> func){
              result_type res(func.name,func.expr);
              handle_function_arg(res,func.t1,"_1_",access_names_map);
              handle_function_arg(res,func.t2,"_2_",access_names_map);
              handle_function_arg(res,func.t3,"_3_",access_names_map);
              handle_function_arg(res,func.t4,"_4_",access_names_map);
              handle_function_arg(res,func.t5,"_5_",access_names_map);
              return res;


          }
      };

      template<class T>
      typename dummy2exptree_impl<T>::result_type dummy2exptree(std::map<viennacl::backend::mem_handle,shared_infos> & access_names_map,
                                                                T const & t){
          return dummy2exptree_impl<T>::execute(access_names_map,t);
      }




  } // namespace generator
} // namespace viennacl


#endif
