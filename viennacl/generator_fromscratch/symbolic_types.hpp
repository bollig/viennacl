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
          typedef typename LHS::ScalarType ScalarType;
          viennacl::backend::mem_handle const & handle() const{
              return handle_;
          }

      public:
          inprod_infos(shared_infos_map_t & shared_infos,
                       temporaries_map_t & temporaries_map,
                       LHS const & lhs, RHS const & rhs) :
              inprod_infos_base(new LHS(lhs), new RHS(rhs)
                                ,new step_t(inprod_infos_base::compute)){
              temporaries_map_t::iterator it = temporaries_map.insert(std::make_pair(this,viennacl::backend::mem_handle())).first;
              viennacl::backend::memory_create(it->second,sizeof(ScalarType)*128);
              handle_ = it->second;
              infos_ = &shared_infos.insert(std::make_pair(it->second,shared_infos_t(shared_infos.size(),print_type<ScalarType>::value()))).first->second;
          }

          std::string kernel_arguments() const{
              return "__global " + scalartype() + "*" + " " + name();
          }
      private:
          viennacl::backend::mem_handle handle_;
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

        kernel_argument& get(){
            static self_type res;
            return res;
        }

    };

      /**
      * @brief Symbolic vector type
      *
      * @tparam SCALARTYPE The Scalartype of the vector in the generated code
      * @tparam ALIGNMENT The Alignment of the vector in the generated code
      */
      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      class symbolic_vector : public vec_infos_base{
        private:
          typedef symbolic_vector<SCALARTYPE,ALIGNMENT> self_type;
        public:
          typedef viennacl::vector<SCALARTYPE,ALIGNMENT> vcl_vec_t;
          typedef SCALARTYPE ScalarType;
          symbolic_vector(shared_infos_map_t & map
                          ,vcl_vec_t const & vcl_vec) : vcl_vec_(vcl_vec){
            infos_= &map.insert(std::make_pair(vcl_vec_.handle(),shared_infos_t(map.size(),print_type<ScalarType>::value()))).first->second;
          }
          virtual viennacl::backend::mem_handle const & handle() const{ return vcl_vec_.handle(); }

        private:
          vcl_vec_t const & vcl_vec_;
      };


      /**
      * @brief Symbolic matrix type
      *
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

      template<class LHS1, class RHS1, class LHS2, class RHS2>
      struct get_symbolic_type<inner_prod_wrapper<LHS1, RHS1>, LHS2, RHS2> { typedef inprod_infos<LHS2,RHS2> type; };
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

          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     T const & t){
              return result_type(dummy2exptree_impl<Lhs>::execute(shared_infos, temporaries, t.lhs())
                                 ,dummy2exptree_impl<Rhs>::execute(shared_infos, temporaries, t.rhs()));
          }
      };

      template<class ScalarType, unsigned int Alignment>
      struct dummy2exptree_impl<dummy_vector<ScalarType, Alignment> >{
          typedef symbolic_vector<ScalarType, Alignment> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries_,
                                     dummy_vector<ScalarType,Alignment> const & v){
              return result_type(shared_infos, v.vec());
          }
      };

      template<class LHS, class RHS>
      struct dummy2exptree_impl<inner_prod_wrapper<LHS,RHS> >{
      private:
          typedef typename dummy2exptree_impl<LHS>::result_type LhsResult;
          typedef typename dummy2exptree_impl<RHS>::result_type RhsResult;
      public:
          typedef inprod_infos<LhsResult,RhsResult> result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     inner_prod_wrapper<LHS,RHS> const & v){
              return result_type(shared_infos, temporaries,
                                 dummy2exptree_impl<LHS>::execute(shared_infos,temporaries,v.lhs()),
                                 dummy2exptree_impl<RHS>::execute(shared_infos,temporaries,v.rhs()));
          }
      };




      template<class T1, class T2, class T3, class T4, class T5>
      struct dummy2exptree_impl<function_wrapper_impl<T1,T2,T3,T4,T5> >{
      private:
          template<class U>
          static void handle_function_arg(symbolic_function & fun, U const * t, std::string name
                              , shared_infos_map_t & shared_infos
                              , temporaries_map_t & temporaries)
          {

              fun.add_arg(name, dummy2exptree_impl<U>::execute(shared_infos,temporaries,*t));
          }

          static void handle_function_arg(symbolic_function & fun, void const* t, std::string name
                              , shared_infos_map_t & shared_infos
                              , temporaries_map_t & temporaries)
          { }

      public:
          typedef symbolic_function result_type;
          static result_type execute(shared_infos_map_t & shared_infos,
                                     temporaries_map_t & temporaries,
                                     function_wrapper_impl<T1,T2,T3,T4,T5> func){
              result_type res(func.name,func.expr);
              handle_function_arg(res,func.t1,"_1_",shared_infos,temporaries);
              handle_function_arg(res,func.t2,"_2_",shared_infos,temporaries);
              handle_function_arg(res,func.t3,"_3_",shared_infos,temporaries);
              handle_function_arg(res,func.t4,"_4_",shared_infos,temporaries);
              handle_function_arg(res,func.t5,"_5_",shared_infos,temporaries);
              return res;


          }
      };

      template<class T>
      typename dummy2exptree_impl<T>::result_type dummy2exptree(shared_infos_map_t & shared_infos,
                                                                temporaries_map_t & temporaries,
                                                                T const & t){
          return dummy2exptree_impl<T>::execute(shared_infos,temporaries,t);
      }




  } // namespace generator
} // namespace viennacl


#endif
