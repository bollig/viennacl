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
          template<class SharedInfosMapT, class TemporariesMapT>
          inprod_infos(SharedInfosMapT & shared_infos,
                       TemporariesMapT & temporaries_map,
                       LHS const & lhs, RHS const & rhs) :
              inprod_infos_base(new LHS(lhs), new RHS(rhs)
                                ,new step_t(inprod_infos_base::compute)){
              typename TemporariesMapT::iterator it = temporaries_map.insert(std::make_pair(this,viennacl::backend::mem_handle())).first;
              viennacl::backend::memory_create(it->second,sizeof(ScalarType)*128);
              handle_ = it->second;
              infos_ = &shared_infos.insert(std::make_pair(it->second,shared_infos_t(shared_infos.size(),print_type<ScalarType>::value()))).first->second;
          }

           void enqueue(unsigned int & arg, viennacl::ocl::kernel & k) const{
               k.arg(arg++,handle_.opencl_handle());
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
          template<class SharedInfosMapT>
          symbolic_vector(SharedInfosMapT & map
                          ,vcl_vec_t const & vcl_vec) : vcl_vec_(vcl_vec){
            infos_= &map.insert(std::make_pair(vcl_vec_.handle(),shared_infos_t(map.size(),print_type<ScalarType>::value()))).first->second;
          }
          virtual viennacl::backend::mem_handle const & handle() const{ return vcl_vec_.handle(); }
          void enqueue(unsigned int & n_arg, viennacl::ocl::kernel & k) const{
              k.arg(n_arg++,vcl_vec_);
              k.arg(n_arg++,cl_uint(vcl_vec_.internal_size()));
          }

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





  } // namespace generator
} // namespace viennacl


#endif
