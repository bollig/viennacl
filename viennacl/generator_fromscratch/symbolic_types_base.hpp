#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_BASE_HPP


#include "viennacl/forwards.h"
#include "viennacl/ocl/utils.hpp"

#include "viennacl/generator_fromscratch/utils.hpp"
#include "viennacl/generator_fromscratch/forwards.h"

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
            shared_infos_t(unsigned int id, std::string scalartype, unsigned int alignment = 1) : id_(id), name_("arg"+to_string(id)), scalartype_(scalartype), alignment_(alignment){ }
            std::string const & access_name(unsigned int i){ return access_names_.at(i); }
            void  access_name(unsigned int i, std::string const & name_){ access_names_[i] = name_; }
            std::string const & name() const{ return name_; }
            unsigned int id() const{ return id_; }
            std::string const & scalartype() const{ return scalartype_; }
            unsigned int alignment() const{ return alignment_; }
            void alignment(unsigned int val) { alignment_ = val; }
        private:
            std::map<unsigned int,std::string> access_names_;
            unsigned int id_;
            std::string name_;
            std::string scalartype_;
            unsigned int alignment_;
        };







        class kernel_argument;

        class infos_base{
        public:
            typedef std::pair<unsigned long, unsigned int> id_res_t;
            virtual std::string generate(unsigned int i) const { return ""; }
            virtual id_res_t id() const = 0;
            virtual ~infos_base(){ }
        };


        class op_infos_base : public infos_base{
        public:
            id_res_t id() const { return std::make_pair(id_,4); }
            bool is_assignment() const { return is_assignment_; }
        protected:
            op_infos_base(unsigned long id, bool is_assignment) : id_(id), is_assignment_(is_assignment){ }
        private:
            unsigned long id_;
            bool is_assignment_;
        };

        class nonarithmetic_op_infos_base : public op_infos_base{
        public:
            nonarithmetic_op_infos_base(unsigned int id) : op_infos_base(id,false){ }
        };

        class arithmetic_op_infos_base : public op_infos_base{
        public:
            std::string generate(unsigned int i) const{ return expr_; }
        protected:
            arithmetic_op_infos_base( std::string const & expr, std::string const & name, bool is_assignment, unsigned long id) :  op_infos_base(id,is_assignment), expr_(expr), name_(name){ }
        private:
            std::string expr_;
            std::string name_;
        };

        class assign_type : public arithmetic_op_infos_base{
        public:
            assign_type() : arithmetic_op_infos_base(" = ", "eq",true,0){ }
        };

        class add_type : public arithmetic_op_infos_base{
        public:
          add_type() : arithmetic_op_infos_base(" + ", "p",false,1){ }
        };

        class inplace_add_type : public arithmetic_op_infos_base{
        public:
          inplace_add_type() : arithmetic_op_infos_base(" += ", "p_eq",true,2){ }
        };

        class sub_type : public arithmetic_op_infos_base{
        public:
          sub_type() : arithmetic_op_infos_base(" - ", "m",false,3){ }
        };

        class inplace_sub_type : public arithmetic_op_infos_base{
        public:
          inplace_sub_type() : arithmetic_op_infos_base(" -= ", "m_eq",true,4){ }
        };

        class scal_mul_type : public arithmetic_op_infos_base{
        public:
          scal_mul_type() : arithmetic_op_infos_base(" * ", "mu",false,5){ }
        };

        class inplace_scal_mul_type : public arithmetic_op_infos_base{
            inplace_scal_mul_type() : arithmetic_op_infos_base(" *= ", "mu_eq",true,6){ }
        };


        class scal_div_type : public arithmetic_op_infos_base{
          scal_div_type() : arithmetic_op_infos_base(" / ", "div", false,7){ }
        };

        class inplace_scal_div_type :  public arithmetic_op_infos_base{
            inplace_scal_div_type() : arithmetic_op_infos_base(" /= ", "div_eq", true,8){ }
        };

        class inner_prod_type : public nonarithmetic_op_infos_base{
        public:
            inner_prod_type() : nonarithmetic_op_infos_base(9){ }
        };

        class binary_tree_infos_base{
        public:
            infos_base & lhs() const{ return *lhs_; }
            infos_base & rhs() const{ return *rhs_; }
            op_infos_base & op() { return *op_; }
            infos_base::id_res_t id() const {
               infos_base::id_res_t lhs_res = lhs_->id();
               infos_base::id_res_t op_res = op_->id();
               infos_base::id_res_t rhs_res = rhs_->id();
               return std::make_pair(lhs_res.first | op_res.first<<lhs_res.second | rhs_res.first << (lhs_res.second + op_res.second)
                                                   ,lhs_res.second+op_res.second+rhs_res.second);
            }
        protected:
            binary_tree_infos_base(infos_base * lhs, op_infos_base * op, infos_base * rhs) : lhs_(lhs), op_(op), rhs_(rhs){        }
            viennacl::tools::shared_ptr<infos_base> lhs_;
            viennacl::tools::shared_ptr<op_infos_base> op_;
            viennacl::tools::shared_ptr<infos_base> rhs_;
        };


        class unary_tree_infos_base{
        public:
            infos_base & sub(){ return *sub_; }

            std::string generate(unsigned int i) const{
                return "(-" + sub_->generate(i) +")";
            }

        protected:
            unary_tree_infos_base(infos_base * sub) : sub_(sub){        }
            viennacl::tools::shared_ptr<infos_base> sub_;
        };


        class arithmetic_tree_infos_base :  public infos_base,public binary_tree_infos_base{
        public:
            std::string generate(unsigned int i) const { return "(" + lhs_->generate(i) + op_->generate(i) + rhs_->generate(i) + ")"; }
            arithmetic_tree_infos_base( infos_base * lhs, arithmetic_op_infos_base* op, infos_base * rhs) :  binary_tree_infos_base(lhs,op,rhs){        }
            virtual id_res_t id() const { return binary_tree_infos_base::id(); }
        private:
        };

        class vector_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            vector_expression_infos_base( infos_base * lhs, arithmetic_op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class scalar_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            scalar_expression_infos_base( infos_base * lhs, arithmetic_op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class matrix_expression_infos_base : public arithmetic_tree_infos_base{
        public:
            matrix_expression_infos_base( infos_base * lhs, arithmetic_op_infos_base* op, infos_base * rhs) : arithmetic_tree_infos_base( lhs,op,rhs){ }
        };

        class kernel_argument : public infos_base{
        public:
            virtual viennacl::backend::mem_handle const & handle() const = 0;
            kernel_argument( ) { }
            void access_name(unsigned int i, std::string const & new_name) { infos_->access_name(i,new_name); }
            virtual ~kernel_argument(){ }
            virtual std::string generate(unsigned int i) const { return infos_->access_name(i); }
            std::string name() const { return infos_->name(); }
            std::string const & scalartype() const { return infos_->scalartype(); }
            std::string aligned_scalartype() const {
                unsigned int alignment = infos_->alignment();
                std::string const & scalartype = infos_->scalartype();
                if(alignment==1){
                    return scalartype;
                }
                else{
                    assert( (alignment==2 || alignment==4 || alignment==8 || alignment==16) && "Invalid alignment");
                    return scalartype + to_string(alignment);
                }
            }
            unsigned int alignment() const { return infos_->alignment(); }
            void alignment(unsigned int val) { infos_->alignment(val); }
            virtual std::string arguments_string() const = 0;
            virtual void enqueue(unsigned int & arg, viennacl::ocl::kernel & k) const = 0;
        protected:
            shared_infos_t* infos_;
        };

        class user_kernel_argument : public kernel_argument{
        public:
            user_kernel_argument(typename id_res_t::second_type id) : id_(id){ }
            id_res_t id() const { return std::make_pair(6,id_); }
        private:
            typename id_res_t::second_type id_;
        };

        class temporary_kernel_argument : public kernel_argument{ };

        class cpu_scal_infos_base : public user_kernel_argument{
        protected:
            cpu_scal_infos_base(unsigned int id) : user_kernel_argument(id) { }
        public:
            virtual std::string arguments_string() const{
                return scalartype() + " " + name();
            }
        };

        class gpu_scal_infos_base : public user_kernel_argument{
        protected:
            gpu_scal_infos_base(unsigned int id) : user_kernel_argument(id) { }
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
            std::string generate(unsigned int i) const{
                if(*step_==compute){
                    return sum_name() + " += " "dot((" + lhs_->generate(i) +  ")" " , " "(" + rhs_->generate(i) + "))" ;
                }
                return infos_->access_name(0);
            }
            std::string sum_name() const{
                return name()+"_sum";
            }

            unsigned int n_groupsize_used_for_compute(){ return 0; }

        protected:
            inprod_infos_base(infos_base * lhs
                              , infos_base * rhs
                              ,step_t * step): binary_tree_infos_base(lhs,new inner_prod_type(),rhs),step_(step){

            }


        private:
            viennacl::tools::shared_ptr<step_t> step_;
        };

        class vec_infos_base : public user_kernel_argument{
        public:
            std::string  size() const{ return name() + "_size"; }
//            std::string  internal_size() const{ return name() + "_internal_size";}
            std::string  start() const{ return name() + "_start";}
            std::string  inc() const{ return name() + "_inc";}
            std::string arguments_string() const{ return  " __global " + aligned_scalartype() + "*"  + " " + name()
                                                                     + ", unsigned int " + size();                                                                     ;
                                                }
            virtual ~vec_infos_base(){ }
        protected:
            vec_infos_base(unsigned int id) : user_kernel_argument(id) { }
        };


        class mat_infos_base : public user_kernel_argument{
        public:
            std::string  size1() const{ return name() +"size1_"; }
            std::string  size2() const{ return name() +"size2_"; }
            std::string  row_inc() const{ return name() +"row_inc_"; }
            std::string  col_inc() const{ return name() +"col_inc_";}
            std::string  row_start() const{ return name() +"row_start_";}
            std::string  col_start() const{ return name() +"col_start_";}
            std::string arguments_string() const{
                return " __global " + aligned_scalartype() + "*"  + " " + name()
                                                            + ", unsigned int " + row_start()
                                                            + ", unsigned int " + col_start()
                                                            + ", unsigned int " + row_inc()
                                                            + ", unsigned int " + col_inc()
                                                            + ", unsigned int " + size1()
                                                            + ", unsigned int " + size2();
            }
            bool const is_rowmajor() const { return is_rowmajor_; }
            bool const is_transposed() const { return is_transposed_; }
            virtual ~mat_infos_base() { }
        protected:
            mat_infos_base(unsigned int id
                           ,bool is_rowmajor
                           ,bool is_transposed) : user_kernel_argument(id)
                                                  ,is_rowmajor_(is_rowmajor)
                                                  ,is_transposed_(is_transposed){ }
        private:
            bool is_rowmajor_;
            bool is_transposed_;
        };

        class function_base : public infos_base{
        protected:
            typedef std::map<std::string,viennacl::tools::shared_ptr<infos_base> > args_map_t;
        public:
            function_base(std::string const & name) : name_(name), id_(2){ }
            virtual std::string name() const {
                return name_;
            }

            id_res_t id() const{
                typename id_res_t::first_type argres1=0;
                typename id_res_t::second_type argres2=0;
                for(args_map_t::const_iterator it = args_map_.begin() ; it != args_map_.end() ; ++it){
                    id_res_t current_res = it->second->id();
                    argres1 |=  current_res.first << argres2;
                    argres2 += current_res.second;
                }
                return std::make_pair(0 | id_<<6 | argres1<<10
                                      ,6 + 10 + argres2);
            }

            std::list<infos_base*> args() const{
                std::list<infos_base*> res;
                for(args_map_t::const_iterator it = args_map_.begin() ; it!= args_map_.end() ; ++it)
                    res.push_back(it->second.get());
                return res;
            }

        protected:
            std::string name_;
            args_map_t args_map_;
            unsigned int id_;
        };

        template<class T1, class T2=void, class T3=void>
        class symbolic_function : public function_base{
        public:
            typedef typename T1::ScalarType ScalarType;

            symbolic_function(std::string const & name,std::string const & expr) : function_base(name), expr_(expr){
            }


            template<class T>
            void add_arg(std::string const & arg_name, T const & t){
                args_map_.insert(std::make_pair(arg_name, new T(t)));
            }


            virtual std::string generate(unsigned int i) const {
                std::string res(expr_);
                for(args_map_t::const_iterator it = args_map_.begin() ; it!= args_map_.end() ; ++it)
                    replace_all_occurences(res,it->first,it->second->generate(i));
                return res;
            }


        private:
            std::string expr_;
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


//        std::pair<unsigned long, unsigned int> get_id(infos_base const * arg){
//            typedef std::pair<unsigned long, unsigned int> res_t;
//            if(op_infos_base* op_ptr = dynamic_cast<op_infos_base*>(&arg)){
//                if
//            }
//            else if(user_kernel_argument* user_kernel_arg_ptr = dynamic_cast<user_kernel_argument*>(&arg)){

//            }
//            else if(arithmetic_tree_infos_base* ar_tree_ptr = dynamic_cast<arithmetic_tree_infos_base*>(&arg)){
//                res_t lhs_res = get_id(*ar_tree_ptr->lhs());
//                res_t op_res = get_id(*ar_tree_ptr->op());
//                res_t rhs_res = get_id(*ar_tree_ptr->rhs());
//                return std::make_pair(lhs_res.first | op_res.first<<lhs_res.second | rhs_res.first << (lhs_res.second + op_res.second)
//                                       ,lhs_res.second+op_res.second+rhs_res.second);
//            }
//            else if(function_base* func_ptr = dynamic_cast<function_base*>(&arg)){

//            }
//            return 0;
//        }

        template<class T, class Pred>
        void extract_as(infos_base* root, std::set<T*, deref_less> & args, Pred pred){
            if(arithmetic_tree_infos_base* p = dynamic_cast<arithmetic_tree_infos_base*>(root)){
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
