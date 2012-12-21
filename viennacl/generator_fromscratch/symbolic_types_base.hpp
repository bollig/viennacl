#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP

#include "viennacl/generator_fromscratch/utils.hpp"
#include "viennacl/generator_fromscratch/forwards.h"
#include "viennacl/backend/memory.hpp"
#include <map>
#include <set>
#include <list>

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

        struct shared_infos_t{
        public:
            shared_infos_t(unsigned int id, std::string scalartype) : id_(id), name_("arg"+to_string(id)), scalartype_(scalartype){ }
            std::string & access_name(){ return access_name_; }
            std::string const & name() const{ return name_; }
            unsigned int id() const{ return id_; }
            std::string const & scalartype() const{ return scalartype_; }
        private:
            std::string access_name_;
            unsigned int id_;
            std::string name_;
            std::string scalartype_;
        };



        class op_infos_base{
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
            infos_base & lhs() const{ return *lhs_; }
            infos_base & rhs() const{ return *rhs_; }


        protected:
            binary_tree_infos_base(infos_base * lhs, infos_base * rhs) : lhs_(lhs), rhs_(rhs){        }
            viennacl::tools::shared_ptr<infos_base> lhs_;
            viennacl::tools::shared_ptr<infos_base> rhs_;
        };

        class kernel_argument;

        class infos_base{
        public:
            virtual std::string generate() const = 0;
            virtual std::string name() const = 0;
            virtual ~infos_base(){ }
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
        public:
            virtual viennacl::backend::mem_handle const & handle() const = 0;
            kernel_argument( ) { }
            void access_name(std::string const & new_name) { infos_->access_name() = new_name; }
            virtual ~kernel_argument(){ }
            virtual std::string generate() const { return infos_->access_name(); }
            virtual std::string name() const { return infos_->name(); }
            std::string const & scalartype() const { return infos_->scalartype(); }
            virtual std::string arguments_string() const = 0;
        protected:
            shared_infos_t* infos_;
        };

        class user_kernel_argument : public kernel_argument{ };

        class temporary_kernel_argument : public kernel_argument{ };

        class cpu_scal_infos_base : public user_kernel_argument{
        protected:
            cpu_scal_infos_base() { }
        public:
            virtual std::string arguments_string() const{
                return scalartype() + " " + name();
            }
        };

        class gpu_scal_infos_base : public user_kernel_argument{
        protected:
            gpu_scal_infos_base() { }
        public:
            virtual std::string arguments_string() const{
                return  "__global " + scalartype() + "*"  + " " + name();
            }
        };

        class inprod_infos_base : public binary_tree_infos_base,public temporary_kernel_argument{
        public:
            enum step_t{compute,reduce};
            step_t step(){ return *step_; }
            void step(step_t s){ *step_ = s; }
            local_memory make_local_memory(unsigned int size){
                return local_memory(name()+"_local",size,scalartype());
            }
            std::string arguments_string() const{
                return "__global " + scalartype() + "*" + " " + name();
            }
            std::string generate() const{
                if(*step_==compute){
                    return sum_name() + " += " "(" + lhs_->generate() +  ")" " * " "(" + rhs_->generate() + ")" ;
                }
                return infos_->access_name();
            }
            std::string sum_name() const{
                return name()+"_sum";
            }

        protected:
            inprod_infos_base(infos_base * lhs
                              , infos_base * rhs
                              ,step_t * step): binary_tree_infos_base(lhs,rhs),step_(step){

            }


        private:
            viennacl::tools::shared_ptr<step_t> step_;
        };

        class vec_infos_base : public user_kernel_argument{
        public:
            std::string  size() const{ return name() + "_size_"; }
            std::string  internal_size() const{ return name() + "_internal_size_";}
            std::string  start() const{ return name() + "_start";}
            std::string  inc() const{ return name() + "_inc";}
            std::string arguments_string() const{ return  " __global " + scalartype() + "*"  + " " + name()
                                                                     + ", unsigned int " + size()
                                                                     + ", unsigned int " + internal_size();}
            virtual ~vec_infos_base(){ }
        protected:
            vec_infos_base() { }
        };


        class mat_infos_base : public user_kernel_argument{
        public:
            std::string  size1() const{ return name() +"size1_"; }
            std::string  size2() const{ return name() +"size2_"; }
            std::string  internal_size1() const{ return name() +"internal_size1_"; }
            std::string  internal_size2() const{ return name() +"internal_size2_"; }
            std::string  row_inc() const{ return name() +"row_inc_"; }
            std::string  col_inc() const{ return name() +"col_inc_";}
            std::string  row_start() const{ return name() +"row_start_";}
            std::string  col_start() const{ return name() +"col_start_";}
            std::string arguments_string() const{
                " __global " + scalartype() + "*"  + " " + name()
                                                            + ", unsigned int " + row_start()
                                                            + ", unsigned int " + col_start()
                                                            + ", unsigned int " + row_inc()
                                                            + ", unsigned int " + col_inc()
                                                            + ", unsigned int " + size1()
                                                            + ", unsigned int " + size2()
                                                            + ", unsigned int " + internal_size1()
                                                            + ", unsigned int " + internal_size2();
            }
            bool const is_rowmajor() const { return is_rowmajor_; }
            bool const is_transposed() const { return is_transposed_; }
            virtual ~mat_infos_base() { }
        protected:
            mat_infos_base(bool is_rowmajor
                           ,bool is_transposed) : is_rowmajor_(is_rowmajor)
                                                  ,is_transposed_(is_transposed){ }
        private:
            bool is_rowmajor_;
            bool is_transposed_;
        };

        class function_base : public infos_base{
        protected:
            typedef std::map<std::string,viennacl::tools::shared_ptr<infos_base> > args_map_t;
        public:
            function_base(std::string const & name) : name_(name){ }
            virtual std::string name() const {
                return name_;
            }

            virtual std::list<infos_base*> args() const = 0;

        protected:
            std::string name_;

        };

        class symbolic_function : public function_base{
        public:
            symbolic_function(std::string const & name,std::string const & expr) : function_base(name), expr_(expr){
            }


            template<class T>
            void add_arg(std::string const & arg_name, T const & t){
                args_map_.insert(std::make_pair(arg_name, new T(t)));
            }


            std::list<infos_base*> args() const{
                std::list<infos_base*> res;
                for(args_map_t::const_iterator it = args_map_.begin() ; it!= args_map_.end() ; ++it)
                    res.push_back(it->second.get());
                return res;
            }

            virtual std::string generate() const {
                std::string res(expr_);
                for(args_map_t::const_iterator it = args_map_.begin() ; it!= args_map_.end() ; ++it)
                    replace_all_occurences(res,it->first,it->second->generate());
                return res;
            }


        private:
            std::string expr_;
            args_map_t args_map_;
        };


        bool operator<(infos_base const & first, infos_base const & other){
            if(binary_tree_infos_base const * t = dynamic_cast<binary_tree_infos_base const *>(&first)){
                return t->lhs() < other || t->rhs() < other;
            }
            else if(binary_tree_infos_base const * p= dynamic_cast<binary_tree_infos_base const *>(&other)){
                return first < p->lhs() || first < p->rhs();
            }
            else if(user_kernel_argument const * t = dynamic_cast<user_kernel_argument const *>(&first)){
                  if(user_kernel_argument const * p = dynamic_cast<user_kernel_argument const*>(&other)){
                     return t->handle() < p->handle();
                  }
            }
            return false;
       }

        typedef std::map<viennacl::backend::mem_handle, shared_infos_t> shared_infos_map_t;
        typedef std::map<kernel_argument*,viennacl::backend::mem_handle,deref_less> temporaries_map_t;

        template<class T, class Pred>
        void extract_as(infos_base* root, std::set<T*, deref_less> & args, Pred pred){
            if(user_kernel_argument* p = dynamic_cast<user_kernel_argument*>(root)){

            }
            else if(arithmetic_tree_infos_base* p = dynamic_cast<arithmetic_tree_infos_base*>(root)){
                extract_as(&p->lhs(), args,pred);
                extract_as(&p->rhs(),args,pred);
            }
            else if(function_base* p = dynamic_cast<function_base*>(root)){
                std::list<infos_base*> func_args(p->args());
                for(std::list<infos_base*>::const_iterator it = func_args.begin(); it!= func_args.end(); ++it){
                    extract_as(*it,args,pred);
                }
            }
            else if(inprod_infos_base* p = dynamic_cast<inprod_infos_base*>(root)){
                if(p->step() == inprod_infos_base::compute){
                    extract_as(&p->lhs(), args,pred);
                    extract_as(&p->rhs(),args,pred);
                }
            }
            if(T* t = dynamic_cast<T*>(root))
                if(pred(t)) args.insert(t);
        }

    }

}
#endif // SYMBOLIC_TYPES_BASE_HPP
