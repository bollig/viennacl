#ifndef VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP

#include "viennacl/generator_fromscratch/symbolic_types_base.hpp"
#include "viennacl/generator_fromscratch/code_generation/utils.hpp"
#include <algorithm>

namespace viennacl{

    namespace generator{

        namespace code_generation{

                struct lid_has_same_size : public std::binary_function<local_memory,local_memory,bool>{
                    bool operator()(local_memory const & lhs, local_memory const & rhs) const{ return lhs.size() == rhs.size(); }
                };


                class optimization_profile{
                protected:
                    typedef unsigned int size_type;

                    virtual void print(std::ostream & os) const{
                        os << " Alignment : " << alignment_ << " | "
                           << " Unroll : " << loop_unroll_ << " | "
                           << " Local Work Size 0 : " << local_work_size_[0] << " | "
                           << " Global Work Size 0 : " << global_work_size_[0] << " | "
                           << " Local Work Size 1 : " << local_work_size_[1] << " | "
                           << " Global Work Size 1 : " << global_work_size_[1] ;
                    }

                public:

                    optimization_profile() : alignment_(1), loop_unroll_(1){
                        local_work_size_[0] = 128;
                        global_work_size_[0] = 128*128;
                    }

                    void load(viennacl::ocl::device const & d){

                    }

                    void config_nd_range(viennacl::ocl::kernel & k) const{
                        k.local_work_size(0,local_work_size_[0]);
                        k.local_work_size(1,local_work_size_[1]);
                        k.global_work_size(0,global_work_size_[0]);
                        k.global_work_size(1,global_work_size_[1]);
                    }

                    void apply(std::list<infos_base*> & expressions){
                        std::set<kernel_argument*,viennacl::generator::deref_less> kernel_arguments;
                        for(std::list<infos_base*>::iterator it = expressions.begin() ; it!=expressions.end() ; ++it){
                            extract_as(*it,kernel_arguments,utils::is_type<kernel_argument>());
                        }
                        for(std::set<kernel_argument*,viennacl::generator::deref_less>::iterator it=kernel_arguments.begin() ; it!= kernel_arguments.end() ; ++it){
                            (*it)->alignment(alignment_);
                        }
                    }

                    unsigned int alignment() const{ return alignment_; }
                    unsigned int loop_unroll(){ return loop_unroll_; }

                    friend std::ostream& operator<<(std::ostream& os, optimization_profile const & );

                    void local_work_size(unsigned int index, size_type val){
                        assert(index==0 || index==1);
                        local_work_size_[index]=val;
                    }

                    size_type local_work_size(unsigned int index) const{
                        assert(index==0 || index==1);
                        return local_work_size_[index];
                    }

                    void global_work_size(unsigned int index, size_type val){
                        assert(index==0 || index==1);
                        global_work_size_[index]=val;
                    }

                    size_type global_work_size(unsigned int index) const{
                        assert(index==0 || index==1);
                        return global_work_size_[index];
                    }

                    virtual ~optimization_profile(){ }

                protected:
                    size_type local_work_size_[2];
                    size_type global_work_size_[2];
                    unsigned int alignment_;
                    unsigned int loop_unroll_;
                };

                class blas1_optimization_profile : public optimization_profile{

                };

                class blas3_optimization_profile : public optimization_profile{
                private:

                    virtual void print(std::ostream & os) const{
                        os << " ML : " << ml_ << " | "
                           << " KL : " << kl_ << " | "
                           << " NL : " << nl_ << " | "
                           << " MS : " << ms_ << " | "
                           << " KS : " << ks_ << " | "
                           << " NS : " << ns_ << " | "
                           << " Alignment : " << alignment_ << " | ";
                    }


                public:

                    blas3_optimization_profile(unsigned int ml, unsigned int kl, unsigned int nl
                                               , unsigned int ms, unsigned int ks, unsigned int ns
                                               , bool use_LHS_shared, bool use_RHS_shared
                                               , unsigned int alignment) : ml_(ml), kl_(kl), nl_(nl), ms_(ms), ks_(ks), ns_(ns), use_LHS_shared_(use_LHS_shared), use_RHS_shared_(use_RHS_shared){
                        alignment_ = alignment;

                        local_work_size_[0] =ml_/ms_;
                        local_work_size_[1] = nl_/ns_;
                    }

                    void set_global_sizes(size_t mat_size1, size_t mat_size2){
                        global_work_size_[0] = (mat_size1+ml_-1)/ml_*local_work_size_[0];
                        global_work_size_[1] = (mat_size2+nl_-1)/nl_*local_work_size_[1];
                    }

                    unsigned int ml() const{ return ml_ ; }
                    unsigned int kl() const{ return kl_ ; }
                    unsigned int nl() const{ return nl_ ; }
                    unsigned int ms() const{ return ms_ ; }
                    unsigned int ks() const{ return ks_ ; }
                    unsigned int ns() const{ return ns_ ; }
                    bool use_LHS_shared() const{ return use_LHS_shared_; }
                    bool use_RHS_shared() const{ return use_RHS_shared_; }

                private:
                    unsigned int ml_;
                    unsigned int kl_;
                    unsigned int nl_;

                    unsigned int ms_;
                    unsigned int ks_;
                    unsigned int ns_;

                    bool use_LHS_shared_;
                    bool use_RHS_shared_;
                };

                std::ostream& operator<<(std::ostream& os, optimization_profile const & prof){
                    prof.print(os);
                    return os;
                }

                void compute_reductions_samesize(utils::kernel_generation_stream& kss, std::list<local_memory> const & lmems){
                   //Same size assumption
                   assert(std::adjacent_find(lmems.begin(), lmems.end(), std::not2(lid_has_same_size()))==lmems.end() && " Calling the wrong function for reducing inner products of different sizes! ");
                   unsigned int size = lmems.front().size();
                   for(unsigned int stride = size/2 ; stride>0 ; stride /=2){
                       kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                       kss << "if(get_local_id(0) < " << to_string(stride) << "){" << std::endl;
                       kss.inc_tab();
                       for(std::list<local_memory>::const_iterator it = lmems.begin(); it != lmems.end() ; ++it){
                           kss <<  it->access("get_local_id(0)") <<  " += " <<  it->access("get_local_id(0)+" + to_string(stride)) << ";" << std::endl;
                       }
                       kss.dec_tab();
                       kss << "}" << std::endl;
                   }
                }

                class blas1_generator{
                private:
                    struct is_compute :  public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const { return p->step() == inprod_infos_base::compute; } };
                    struct is_reduce : public std::unary_function<inprod_infos_base*, bool>{ bool  operator()(inprod_infos_base* p) const{ return p->step() == inprod_infos_base::reduce; } };

                public:
                    blas1_generator(std::list<infos_base * > const & vector_expressions
                                    , std::list<infos_base * > const & matrix_expressions
                                    , std::list<infos_base * > const & scalar_expressions
                                    , optimization_profile * kernel_config): vector_expressions_(vector_expressions), matrix_expressions_(matrix_expressions), scalar_expressions_(scalar_expressions), optimization_profile_(kernel_config)
                    {
                        std::set<inprod_infos_base*> inprods;
                        for(std::list<infos_base*>::const_iterator it=vector_expressions.begin() ; it!= vector_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                        for(std::list<infos_base*>::const_iterator it=matrix_expressions.begin() ; it!= matrix_expressions.end() ; ++it){
                            extract_as(*it,matrices_,utils::is_type<mat_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                        for(std::list<infos_base*>::const_iterator it=scalar_expressions.begin() ; it!= scalar_expressions.end() ; ++it){
                            extract_as(*it,vectors_,utils::is_type<vec_infos_base>());
                            extract_as(*it,matrices_,utils::is_type<mat_infos_base>());
                            extract_as(*it,gpu_scalars_,utils::is_type<gpu_scal_infos_base>());
                            extract_as(*it,inner_prods_compute_,is_compute());
                            extract_as(*it,inner_prods_reduce_,is_reduce());
                        }
                    }


                    void operator()(utils::kernel_generation_stream& kss){

                        unsigned int local_work_size0 = optimization_profile_->local_work_size(0);
                        unsigned int alignment = optimization_profile_->alignment();
                        unsigned int n_unroll = optimization_profile_->loop_unroll();
                        std::list<vec_infos_base *> assigned_vec;
                        for(std::list<infos_base*>::iterator it=vector_expressions_.begin(); it!= vector_expressions_.end();++it){
                            vector_expression_infos_base* p=static_cast<vector_expression_infos_base*>(*it);
                            if(p->op().is_assignment()==true) assigned_vec.push_back(dynamic_cast<vec_infos_base*>(&p->lhs()));
                        }

                        std::list<mat_infos_base *> assigned_mat;
                        for(std::list<infos_base*>::iterator it=matrix_expressions_.begin(); it!= matrix_expressions_.end();++it){
                            matrix_expression_infos_base* p=static_cast<matrix_expression_infos_base*>(*it);
                            if(p->op().is_assignment()==true) assigned_mat.push_back(dynamic_cast<mat_infos_base*>(&p->lhs()));
                        }

                        std::list<gpu_scal_infos_base *> assigned_scal;
                        for(std::list<infos_base*>::iterator it=scalar_expressions_.begin(); it!= scalar_expressions_.end();++it){
                            if(scalar_expression_infos_base* p=dynamic_cast<scalar_expression_infos_base*>(*it)){
                                if(p->op().is_assignment()==true){
                                    assigned_scal.push_back(dynamic_cast<gpu_scal_infos_base*>(&p->lhs()));
                                }
                            }
                        }

                        code_generation::utils::cache_manager<vec_infos_base> vector_cache(vectors_,assigned_vec,kss);
                        code_generation::utils::cache_manager<mat_infos_base> matrix_cache(matrices_,assigned_mat,kss);
                        code_generation::utils::cache_manager<gpu_scal_infos_base> scalar_cache(gpu_scalars_,assigned_scal,kss);
                        vec_infos_base * first_vector =  NULL;
                        mat_infos_base * first_matrix = NULL;
                        if(vectors_.size())
                            first_vector = *vectors_.begin();
                        if(matrices_.size())
                            first_matrix = *matrices_.begin();
//                        //Assumes same size...

                        if(inner_prods_reduce_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss <<  (*it)->scalartype() <<  " " << (*it)->sum_name()<< " = 0 ;" << std::endl;
                            }
                            kss << "for(unsigned int i = get_global_id(0) ; i <" << optimization_profile_->global_work_size(0)/optimization_profile_->local_work_size(0) << " ; i += get_global_size(0)){" << std::endl;
                            kss.inc_tab();
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                kss << (*it)->sum_name() << " += " << (*it)->name() << "[i];" << std::endl;
                            }
                            kss.dec_tab();
                            kss << "}" << std::endl;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(local_work_size0));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " << (*it)->sum_name()<< ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_reduce_.begin(); it != inner_prods_reduce_.end() ; ++it){
                                (*it)->access_name(0,(*it)->make_local_memory(local_work_size0).access("0"));
                            }
                            kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                        }


                        scalar_cache.fetch_entries(0,"0");
                        if(first_vector){
                            for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                                kss << (*it)->scalartype() << " " << (*it)->sum_name() << " = 0;" << std::endl;
                            }
                            std::list<infos_base*> expressions;
                            std::transform(inner_prods_compute_.begin(), inner_prods_compute_.end(),std::back_inserter(expressions),UnsafeBase2Target<inprod_infos_base,infos_base>());
                            std::copy(vector_expressions_.begin(),vector_expressions_.end(),std::back_inserter(expressions));
                            utils::unroll_gid_loop_contiguous(kss,n_unroll,first_vector->size() + "/" +to_string(alignment) ,expressions,vector_cache);
                        }

                        if(first_matrix){
                            utils::unroll_gid_loop_contiguous(kss,n_unroll,first_matrix->size1()+"*"+first_matrix->size2()+ "/" +to_string(alignment) , matrix_expressions_, matrix_cache);
                        }
                        scalar_cache.writeback_entries(0,"0");

                        if(inner_prods_compute_.size()){
                            std::list<local_memory> local_mems;
                            for( std::set<inprod_infos_base *, viennacl::generator::deref_less>::const_iterator it = inner_prods_compute_.begin(); it != inner_prods_compute_.end() ; ++it){
                                local_mems.push_back((*it)->make_local_memory(local_work_size0));
                                kss << local_mems.back().declare() << ";" << std::endl;
                                kss << local_mems.back().access("get_local_id(0)") << " = " <<  (*it)->sum_name() << ";" << std::endl;
                            }
                            compute_reductions_samesize(kss,local_mems);
                        }
                        for(std::set<inprod_infos_base *, viennacl::generator::deref_less>::iterator it=inner_prods_compute_.begin() ; it!=inner_prods_compute_.end();++it){
                            (*it)->step(inprod_infos_base::reduce);
                            kss << "if(get_local_id(0)==0) " << (*it)->name() << "[get_group_id(0)]" << "=" << (*it)->make_local_memory(local_work_size0).access("0") << ";" << std::endl;
                        }
                    }

                private:
                    std::list<infos_base * >  vector_expressions_;
                    std::list<infos_base * >  matrix_expressions_;
                    std::list<infos_base* > scalar_expressions_;
                    std::set<vec_infos_base *, viennacl::generator::deref_less >  vectors_;
                    std::set<mat_infos_base *, viennacl::generator::deref_less >  matrices_;
                    std::set<gpu_scal_infos_base *, viennacl::generator::deref_less > gpu_scalars_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_compute_;
                    std::set<inprod_infos_base *, viennacl::generator::deref_less > inner_prods_reduce_;
                    optimization_profile * optimization_profile_;


                };

                class blas3_generator{
                    typedef std::set<matmat_prod_infos_base*,viennacl::generator::deref_less> matmat_prods_t;
                public:
                    blas3_generator(std::list<infos_base * > const & blas3_expressions
                                    , blas3_optimization_profile * kernel_config): blas3_expressions_(blas3_expressions), optimization_profile_(kernel_config)
                    {
                        for(std::list<infos_base*>::const_iterator it=blas3_expressions_.begin() ; it!= blas3_expressions_.end() ; ++it){
                            extract_as(*it, matmat_prods_, utils::is_type<matmat_prod_infos_base>());
                            extract_as(*it ,gpu_scalars_,  utils::is_type<gpu_scal_infos_base>());
                            extract_as(*it,matrices_, utils::is_type<mat_infos_base>());
                        }
                    }

                    void operator()(utils::kernel_generation_stream& kss){


                        std::list<mat_infos_base*> assigned;
                        std::set<mat_infos_base*,viennacl::generator::deref_less> lhss;
                        std::set<mat_infos_base*,viennacl::generator::deref_less> rhss;

                        //Fills assigned matrices set
                        for(std::list<infos_base*>::iterator it = blas3_expressions_.begin() ; it!=blas3_expressions_.end(); ++it){
                            matrix_expression_infos_base* p=static_cast<matrix_expression_infos_base*>(*it);
                            if(p->op().is_assignment()) assigned.push_back(static_cast<mat_infos_base*>(&p->lhs()));
                        }

                        //Fills lhs's
                        for(matmat_prods_t::iterator it = matmat_prods_.begin(); it !=matmat_prods_.end(); ++it){
                            extract_as(&(*it)->lhs(),lhss,utils::is_type<mat_infos_base>());
                            extract_as(&(*it)->rhs(),rhss,utils::is_type<mat_infos_base>());
                        }

                        mat_infos_base* first_lhs = *lhss.begin();
                        mat_infos_base* first_rhs = *rhss.begin();
                        mat_infos_base* first_assigned = assigned.front();



                        bool is_lhs_rowmajor = first_lhs->is_rowmajor();
                        bool is_rhs_rowmajor = first_rhs->is_rowmajor();
                        bool is_result_rowmajor = first_assigned->is_rowmajor();

                        bool use_inner_product = is_lhs_rowmajor && !is_rhs_rowmajor;

                        unsigned int alignment = optimization_profile_->alignment();

                        unsigned int ml_res = optimization_profile_->ml();
                        unsigned int nl_res = optimization_profile_->nl();
                        unsigned int ms_res = optimization_profile_->ms();
                        unsigned int ns_res = optimization_profile_->ns();
                        if(is_result_rowmajor){
                            nl_res /= alignment ; ns_res /= alignment;
                        }
                        else{
                            ml_res /= alignment ; ms_res /= alignment;
                        }

                        unsigned int kl_default = optimization_profile_->kl();
                        unsigned int ks_default = optimization_profile_->ks();


                        unsigned int ml_lhs = optimization_profile_->ml();
                        unsigned int kl_lhs = optimization_profile_->kl();
                        unsigned int ms_lhs = optimization_profile_->ms();
                        unsigned int ks_lhs = optimization_profile_->ks();

                        if(first_lhs->is_rowmajor()){
                            kl_lhs /= alignment;
                        }
                        else{
                            ml_lhs /= alignment;
                        }

                        unsigned int kl_rhs = optimization_profile_->kl();
                        unsigned int nl_rhs = optimization_profile_->nl();
                        unsigned int ks_rhs = optimization_profile_->ks();
                        unsigned int ns_rhs = optimization_profile_->ns();
                        if(is_rhs_rowmajor){
                            nl_rhs /= alignment; ns_rhs /= alignment;
                        }
                        else{
                            kl_rhs /= alignment; ks_rhs /= alignment;
                        }

                        bool use_LHS_shared = optimization_profile_->use_LHS_shared();
                        bool use_RHS_shared = optimization_profile_->use_RHS_shared();


                        //Declare local memory
                        if(use_LHS_shared){
                                kss << "__local " << first_lhs->scalartype() << " local_lhs[" << optimization_profile_->ml() * (kl_default+1)<< "];" << std::endl;

                        }

                        std::string res_table_name(first_assigned->name() + "_res");
                        for(unsigned int m=0; m< ms_res; ++m){
                            for(unsigned int n=0; n < ns_res ; ++n){
                                kss << first_assigned->aligned_scalartype() << " " << res_table_name << "_" << m << "_" << n << " = 0 ;" << std::endl;
                            }
                        }

                        std::string offset_m("get_local_id(0)*" + to_string(ms_lhs));
                        std::string offset_n("get_local_id(1)*" + to_string(ns_rhs));
                        //Helper variables
                        if(is_lhs_rowmajor){
                            kss << "unsigned int aligned_size2_lhs = " << first_lhs->internal_size2() << "/" << alignment<< ";" << std::endl;

                        }
                        else{
                            kss << "unsigned int aligned_size1_lhs = " << first_lhs->internal_size1() << "/" << alignment<< ";" << std::endl;
                            kss << "unsigned int aligned_size2_lhs = " << first_lhs->internal_size2() << ";" << std::endl;
                        }

                        if(is_rhs_rowmajor){
                            kss << "unsigned int aligned_size2_rhs = " << first_rhs->internal_size2() << "/" << alignment<< ";" << std::endl;
                        }
                        else{
                            kss << "unsigned int aligned_size1_rhs = " << first_rhs->internal_size1() << "/" << alignment<< ";" << std::endl;

                        }

                        if(is_result_rowmajor){
                            kss << "unsigned int aligned_size2_res = " << first_assigned->internal_size2() << "/" << alignment<< ";" << std::endl;
                        }
                        else{
                            kss << "unsigned int aligned_size1_res = " << first_assigned->internal_size1() << "/" << alignment<< ";" << std::endl;
                        }

                        kss << "unsigned int block_num = aligned_size2_lhs/" << kl_lhs << ";" << std::endl;

                        std::string res_aligned_size1(first_assigned->internal_size1() + "/" + to_string(first_assigned->alignment()));

                        if(is_result_rowmajor)
                            kss << "__global " << first_assigned->aligned_scalartype() << "* res_ptr = " << first_assigned->name() << " + (get_global_id(0)*" << ms_res << ")*aligned_size2_res + get_global_id(1)*" << ns_res << ";" << std::endl;
                        else
                            kss << "__global " << first_assigned->aligned_scalartype() << "* res_ptr = " << first_assigned->name() << " + (get_global_id(0)*" << ms_res << ") + get_global_id(1)*" << res_aligned_size1 << "*" << ns_res << ";" << std::endl;

                        if(is_rhs_rowmajor)
                            kss << "unsigned int offsetRHS = " << offset_n << " +  get_group_id(1)*" << nl_rhs << ";" << std::endl;
                        else
                            kss << "unsigned int offsetRHS = (" << offset_n << "+  get_group_id(1)*" << nl_rhs << ")*aligned_size1_rhs;" << std::endl;

                        if(is_lhs_rowmajor)
                            kss << "unsigned int offsetLHS = get_group_id(0)*" << "aligned_size2_lhs" << "*" << ml_lhs  << ";" << std::endl;
                        else
                            kss << "unsigned int offsetLHS = get_group_id(0)*" << ml_lhs  << ";" << std::endl;

                        kss << "for(unsigned int bl=0 ; bl<block_num ; ++bl){" << std::endl;
                        kss.inc_tab();
                        if(use_LHS_shared){
                            kss << "__local " << first_lhs->scalartype() << " * ptr_lhs; " << std::endl;
                            kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                            kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << ml_lhs << "; i+= get_local_size(0)){" << std::endl;
                            kss.inc_tab();
                            kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << kl_lhs << "; j+= get_local_size(1)){" << std::endl;
                            kss.inc_tab();
                            if(is_lhs_rowmajor){
                                kss << first_lhs->aligned_scalartype()  << " val_lhs = " << first_lhs->name() <<  "[offsetLHS + j  + aligned_size2_lhs*i];" << std::endl;
                                    kss << " ptr_lhs = local_lhs + i*" << kl_default+1 << "+j*" << alignment<<";" <<std::endl;
                                    for(unsigned int a = 0 ; a < alignment ; ++a){
                                        if(alignment>1)
                                            kss << "*ptr_lhs++ =  val_lhs.s" << a << ";" << std::endl;
                                        else
                                            kss << "*ptr_lhs++ =  val_lhs;" << std::endl;
                                    }
                            }
                            else{
                                kss << first_lhs->aligned_scalartype()  << " val_lhs = " << first_lhs->name() <<  "[offsetLHS + j*aligned_size1_lhs  + i];" << std::endl;
                                kss << " ptr_lhs = local_lhs + i*" << alignment * (kl_default+1) << "+ j;" <<std::endl;
                                for(unsigned int a = 0 ; a < alignment ; ++a){
                                    if(alignment>1)
                                        kss << "*ptr_lhs =  val_lhs.s" << a << ";" << std::endl;
                                    else
                                        kss << "*ptr_lhs =  val_lhs;" << std::endl;
                                    kss << "ptr_lhs += " << kl_default + 1 << ";" << std::endl;
                                }
                            }

                            kss.dec_tab();
                            kss << "}" << std::endl;
                            kss.dec_tab();
                            kss << "}" << std::endl;
                            kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                                for(unsigned int m=0; m<ms_lhs; ++m){
                                    kss << "__local " << first_lhs->scalartype() << "* ptr_lhs_" << m << " = local_lhs + (" << offset_m << "+" << m << ")" << "*" << optimization_profile_->kl()  + 1 << ";" << std::endl;
                                }
                        }
                        else{
                            for(unsigned int m=0; m<ms_lhs; ++m){
                                kss << "__global " << first_lhs->scalartype() << "* ptr_lhs_" << m << " = " << first_lhs->name() << " + offsetLHS + (" << offset_m << "+" << m << ")" << "*aligned_size2_lhs" << ";" << std::endl;
                            }
                        }

                        if(is_rhs_rowmajor){
                            for(unsigned int k = 0 ; k < ks_rhs ; ++k){
                                kss << "__global " << first_rhs->aligned_scalartype() << " * rhs_ptr_" << k<< " = " << first_rhs->name() << " + offsetRHS + aligned_size2_rhs*" << k << ";" << std::endl;
                            }
                        }
                        else{
                            for(unsigned int n = 0 ; n < ns_rhs ; ++n){
                                kss << "__global " << first_rhs->aligned_scalartype() << " * rhs_ptr_" << n << " = " << first_rhs->name() << " + offsetRHS + aligned_size1_rhs*" << n << ";" << std::endl;
                            }
                        }

//                        kss << "unsigned int n_subblock = min(" << kl_default/ks_default << ", (int)(" << "aligned_size2_lhs" << " - bl*" << kl_lhs << ")/ " << ks_lhs << ");" << std::endl;

                        kss << " for(unsigned int bs=0 ; bs < " << kl_default/ks_default  << " ; ++bs){" << std::endl;
                        kss.inc_tab();

                        if(is_rhs_rowmajor){
                            for(unsigned int k = 0 ; k < ks_rhs ; ++k){
                                for(unsigned int n=0 ; n < ns_rhs ; ++n){
                                    kss << first_rhs->aligned_scalartype() << " val_rhs_" << k << "_" << n << " = *rhs_ptr_" << k  << "++;" << std::endl;
                                }
                            }
                        }
                        else{
                            for(unsigned int n=0 ; n < ns_rhs ; ++n){
                                for(unsigned int k = 0 ; k < ks_rhs ; ++k){
                                    kss << first_rhs->aligned_scalartype() << " val_rhs_" << k << "_" << n << " = *rhs_ptr_" << n  << "++;" << std::endl;
                                }
                            }
                        }

                       for(unsigned int k = 0 ; k < ks_lhs ; ++k){
                           for(unsigned int m=0 ; m < ms_lhs ; ++m){
                               kss << first_lhs->scalartype() << " " << "val_lhs_" << m << "_" << k << " = " << "* ptr_lhs_" << m << "++;" << std::endl;
                           }
                       }

                        if(is_result_rowmajor){
                            if(use_inner_product){
                                for(unsigned int m=0 ; m < ms_lhs ; ++m){
                                    for(unsigned int n=0 ; n < ns_rhs ; ++n){
                                        for(unsigned int k = 0 ; k < ks_default/alignment ; ++k){
                                            if(alignment>1){
                                                kss << res_table_name<< "_"<<m<<"_" << n/alignment << ".s" << n%alignment<< " += dot(";
                                                kss << "(" << first_lhs->aligned_scalartype() << ")(" ;
                                                for(unsigned int a = 0; a < alignment; ++a){
                                                    kss << "val_lhs_" << m << "_" << k*alignment+a;
                                                    if(a<alignment-1) kss << ",";
                                                }
                                                kss << "), " << "val_rhs_" << k << "_" << n << ");" << std::endl;
                                            }
                                            else
                                                kss << res_table_name<< "_"<<m<<"_" << n << " += dot(val_lhs_" << m << "_" << k << " , " << "val_rhs_" << k << "_" << n << ");" << std::endl;

                                        }
                                    }
                                }
                            }
                            else{
                                for(unsigned int k = 0 ; k < ks_default ; ++k){
                                    for(unsigned int n=0 ; n < ns_rhs ; ++n){
                                        for(unsigned int m=0 ; m < ms_lhs ; ++m){
                                                kss << res_table_name<< "_"<<m<<"_" << n << " += " ;
                                                kss << "val_lhs_" << m << "_" << k;
                                                kss << "*";
                                                kss <<" val_rhs_" << k << "_" << n << ";" << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                        else{
                            for(unsigned int k = 0 ; k < ks_default ; ++k){
                                for(unsigned int n=0 ; n < ns_res ; ++n){
                                    for(unsigned int m=0 ; m < ms_lhs ; ++m){
                                        kss << res_table_name<< "_"<<m/alignment<<"_" << n << ".s" << m%alignment << " += " ;
                                        kss << "val_lhs_" << m << "_" << k;
                                        kss << "*";
                                        kss <<" val_rhs_" << k << "_" << n/alignment << ".s" << n%alignment << ";" << std::endl;
                                    }
                                }
                            }
                        }

                        if(is_rhs_rowmajor){
                            for(unsigned int k=0 ; k<ks_default ; ++k){
                                kss << "rhs_ptr_" << k << " += " << ks_rhs << "*aligned_size2_rhs - " << ns_rhs << ";" << std::endl;
                            }
                        }

                        kss.dec_tab();
                        kss << "}" << std::endl;
                        if(is_rhs_rowmajor){
                            kss << "offsetRHS += aligned_size2_rhs" << "*" << kl_rhs << ";" << std::endl;
                        }
                        else{
                            kss << "offsetRHS += " << kl_rhs << ";" << std::endl;
                        }
                        if(is_lhs_rowmajor)
                            kss << "offsetLHS += " << kl_lhs << ";" << std::endl;
                        else
                            kss << "offsetLHS += " << kl_lhs << "*aligned_size1_lhs;" << std::endl;
                        kss.dec_tab();
                        kss << "}" << std::endl;

                        if(first_assigned->is_rowmajor()){
                            for(unsigned int m=0 ; m < ms_res ; ++m){
                                for(unsigned int n=0 ; n < ns_res ; ++n){
                                    kss << "*res_ptr++=" << res_table_name << "_" << m << "_" << n << ";" << std::endl;
                                }
                                kss << "res_ptr+=aligned_size2_res - " << ns_res << ";" << std::endl;
                            }
                        }
                        else{
                            for(unsigned int n=0 ; n < ns_res ; ++n){
                                for(unsigned int m=0 ; m < ms_res ; ++m){
                                    kss << "*res_ptr++=" << res_table_name << "_" << m << "_" << n << ";" << std::endl;
                                }
                                kss << "res_ptr+=" << res_aligned_size1 << " - " << ms_res << ";" << std::endl;
                            }
                        }
                    }

                private:
                    std::list<infos_base*>  blas3_expressions_;
                    matmat_prods_t matmat_prods_;
                    std::set<mat_infos_base *, viennacl::generator::deref_less >  matrices_;
                    std::set<gpu_scal_infos_base *, viennacl::generator::deref_less > gpu_scalars_;
                    blas3_optimization_profile * optimization_profile_;
                };

        }

    }

}
#endif // BACKEND_HPP
