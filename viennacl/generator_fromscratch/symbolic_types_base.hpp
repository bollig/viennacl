#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP

#include "viennacl/generator_fromscratch/utils.hpp"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/backend/memory.hpp"
namespace viennacl{

    namespace generator{

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
            bool is_assignment() const { return is_assignment_; }
        protected:
            op_infos_base( std::string const & expr, std::string const & name, bool is_assignment) :  expr_(expr), name_(name), is_assignment_(is_assignment){ }
        private:
            std::string expr_;
            std::string name_;
            bool is_assignment_;
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
            infos_base & sub(){ return *sub_; }

            std::string generate() const{
                return "(-" + sub_->generate() +")";
            }

        protected:
            unary_tree_infos_base(infos_base * sub) : sub_(sub){        }
            viennacl::tools::shared_ptr<infos_base> sub_;
        };


        class arithmetic_tree_infos_base :  public infos_base,public binary_tree_infos_base{
        public:
            op_infos_base & op() { return *op_; }
            std::string generate() const { return "(" + lhs_->generate() + op_->generate() + rhs_->generate() + ")"; }
            std::string name() const { return lhs_->name() + op_->name() + rhs_->name(); }
            arithmetic_tree_infos_base( infos_base * lhs, op_infos_base* op, infos_base * rhs) :  binary_tree_infos_base(lhs,rhs), op_(op){        }
        private:
            viennacl::tools::shared_ptr<op_infos_base> op_;
        };

        class vector_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            vector_expression_infos_base( infos_base * lhs, op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class scalar_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            scalar_expression_infos_base( infos_base * lhs, op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class matrix_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            matrix_expression_infos_base( infos_base * lhs, op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class kernel_argument : public infos_base{
        private:
            virtual viennacl::backend::mem_handle handle() const = 0;
        public:
            void access_name(std::string const & new_name) { *access_name_ = new_name; }
            virtual ~kernel_argument(){ }
            virtual std::string generate() const { return *access_name_; }
            virtual std::string name() const { return name_; }
            std::string const & scalartype() const { return scalartype_; }
            kernel_argument( std::string const & name, std::string const & scalartype) :name_(name), scalartype_(scalartype){ }
            virtual std::string kernel_arguments() const = 0;
            bool operator<(kernel_argument const & other) const { return handle() < other.handle(); }
        protected:
            std::string * access_name_;
            std::string name_;
            std::string scalartype_;
        };

        class cpu_scal_infos_base : public kernel_argument{
        protected:
            cpu_scal_infos_base(std::string const & name, std::string const & scalartype) : kernel_argument(name, scalartype){ }
        public:
            std::string kernel_arguments() const{
                return scalartype_ + " " + name_;
            }
        };

        class gpu_scal_infos_base : public kernel_argument{
        protected:
            gpu_scal_infos_base(std::string const & name, std::string const & scalartype) : kernel_argument(name, scalartype){ }
        public:
            std::string kernel_arguments() const{
                return "__global " + scalartype_ + "*"  + " " + name_;
            }
        };

        class inprod_infos_base : public binary_tree_infos_base, public kernel_argument{
        public:
            enum step_t{compute,reduce};
            step_t step(){ return *step_; }
            void step(step_t s){ *step_ = s; }
            local_memory make_local_memory(unsigned int size){
                return local_memory(name_,size,scalartype_);
            }

        protected:
            inprod_infos_base(
                              infos_base * lhs, infos_base * rhs,step_t * step
                              ,std::string const & scalartype): binary_tree_infos_base(lhs,rhs)
                                               ,kernel_argument(lhs->name()+"ip"+rhs->name(),scalartype)
                                              ,step_(step){ }


        private:
            viennacl::tools::shared_ptr<step_t> step_;
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
            vec_infos_base(std::string const & name, std::string const & scalartype) :
                                                            kernel_argument(name,
                                                                            scalartype)
                                                          ,size_(name_+"_size")
                                                          ,internal_size_(""){ }
        private:
            std::string size_;
            std::string internal_size_;
            std::string start_;
            std::string inc_;
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
            mat_infos_base(std::string const & name
                           ,std::string const & scalartype
                           ,bool is_rowmajor
                           ,bool is_transposed) : kernel_argument(name,scalartype)
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


    }

}
#endif // SYMBOLIC_TYPES_BASE_HPP
