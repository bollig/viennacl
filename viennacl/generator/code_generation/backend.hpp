#ifndef VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_BACKEND_HPP

#include "viennacl/generator/symbolic_types_base.hpp"
#include "viennacl/generator/code_generation/utils.hpp"
#include <algorithm>

namespace viennacl{

    namespace generator{

        namespace code_generation{

                struct lid_has_same_size : public std::binary_function<local_memory<1>,local_memory<1>,bool>{
                    bool operator()(local_memory<1> const & lhs, local_memory<1> const & rhs) const{ return lhs.size() == rhs.size(); }
                };


                class optimization_profile{
                protected:
                    typedef unsigned int size_type;


                public:

                    optimization_profile() : alignment_(1), loop_unroll_(1){
                        local_work_size_[0] = 128;
                        global_work_size_[0] = 128*128;
                    }

                    void load(viennacl::ocl::device const & d){

                    }

                    virtual std::string repr() const = 0;

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

//                    friend std::ostream& operator<<(std::ostream& os, optimization_profile const & );

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
                public:
                    virtual std::string repr() const{
                        std::ostringstream oss;
                        oss << " A" << alignment_
                           << " U" << loop_unroll_
                           << " L0" << local_work_size_[0]
                           << " G0" << global_work_size_[0]
                           << " L1" << local_work_size_[1]
                           << " G1" << global_work_size_[1] ;
                        return oss.str();
                    }
                };

                class blas3_optimization_profile : public optimization_profile{
                private:

                public:

                    blas3_optimization_profile(unsigned int ml, unsigned int kl, unsigned int nl
                                               , unsigned int ms, unsigned int ks, unsigned int ns
                                               , bool use_LHS_shared, bool use_RHS_shared
                                               , unsigned int alignment) : ml_(ml), kl_(kl), nl_(nl), ms_(ms), ks_(ks), ns_(ns), use_LHS_shared_(use_LHS_shared), use_RHS_shared_(use_RHS_shared){
                        alignment_ = alignment;

                        local_work_size_[0] =ml_/ms_;
                        local_work_size_[1] = nl_/ns_;
                    }

                    virtual std::string repr() const{
                        std::ostringstream oss;
                        oss << "ML" << ml_
                           << "KL" << kl_
                           << "NL" << nl_
                           << "MS" << ms_
                           << "KS" << ks_
                           << "NS" << ns_
                           << "A" << alignment_
                           << "LHS" << use_LHS_shared_
                           << "RHS" << use_RHS_shared_;
                        return oss.str();
                    }


                    void set_global_sizes(mat_infos_base* mat){
                        local_work_size_[0] =ml_/ms_;
                        global_work_size_[0] = mat->real_size1()/ml_*local_work_size_[0];
                        local_work_size_[1] = nl_/ns_;
                        global_work_size_[1] = mat->real_size2()/nl_*local_work_size_[1];
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

                static std::ostream& operator<<(std::ostream& os, optimization_profile const & prof){
                    os << prof.repr();
                    return os;
                }

                static void compute_reductions_samesize(utils::kernel_generation_stream& kss, std::list<local_memory<1> > const & lmems){
                   //Same size assumption
                   assert(std::adjacent_find(lmems.begin(), lmems.end(), std::not2(lid_has_same_size()))==lmems.end() && " Calling the wrong function for reducing inner products of different sizes! ");
                   unsigned int size = lmems.front().size();
                   for(unsigned int stride = size/2 ; stride>0 ; stride /=2){
                       kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                       kss << "if(get_local_id(0) < " << to_string(stride) << "){" << std::endl;
                       kss.inc_tab();
                       for(std::list<local_memory<1> >::const_iterator it = lmems.begin(); it != lmems.end() ; ++it){
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
                            std::list<local_memory<1> > local_mems;
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
                            std::list<local_memory<1> > local_mems;
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

                    static void transform_block(mat_infos_base const & mat_infos, bool store_shared
                                                , unsigned int & large_block_1, unsigned int & large_block_2
                                                , unsigned int & small_block_1, unsigned int & small_block_2){
                        unsigned int alignment = mat_infos.alignment();
                        if(mat_infos.is_rowmajor()){
                            if(mat_infos.is_transposed()) large_block_1 /= alignment;
                            else large_block_2/=alignment;
                            if(!store_shared){
                                if(mat_infos.is_transposed()) small_block_1/=alignment;
                                else small_block_2/=alignment;
                            }
                        }
                        else{
                            if(mat_infos.is_transposed()) large_block_2 /= alignment;
                            else large_block_1/=alignment;
                            if(!store_shared){
                                if(mat_infos.is_transposed())  small_block_2/=alignment;
                                else    small_block_1/=alignment;
                            }
                        }

                    }


                    static void transform_size(utils::kernel_generation_stream & kss,
                                               mat_infos_base const * mat_infos){
                        if(mat_infos->is_rowmajor()){
                            kss << mat_infos->internal_size2() << "/=" << mat_infos->alignment() << ";" << std::endl;
                        }
                        else{
                            kss << mat_infos->internal_size1() << "/=" << mat_infos->alignment() << ";" << std::endl;
                        }
                    }

                    struct declare_rhs_global_ptr{

                        declare_rhs_global_ptr(utils::kernel_generation_stream & _kss, unsigned int _ks_rhs,unsigned int _ns_rhs,
                                               unsigned int _nl_rhs, std::string const & _offset_n ) : kss(_kss), ks_rhs(_ks_rhs), ns_rhs(_ns_rhs)
                                                                                                        , nl_rhs(_nl_rhs), offset_n(_offset_n)

                        { }

                        void operator()( mat_infos_base * mat) {
                            if(mat->is_rowmajor())
                                for(unsigned int k = 0 ; k < ks_rhs ; ++k){
                                    std::string ptr_name = mat->name() + "_ptr_" + to_string(k);
                                    kss << "__global " << mat->aligned_scalartype() << " * " << ptr_name << " = " << mat->name() << " + " ;
                                    if(mat->is_transposed()) kss<< mat->offset(to_string(k) + " + " + offset_n + " +  get_group_id(1)*" + to_string(nl_rhs),"0");
                                    else kss << mat->offset(to_string(k),offset_n + " +  get_group_id(1)*" + to_string(nl_rhs));
                                    kss << ";" << std::endl;
                                    mat->access_name(k,"*" + ptr_name);
                                }
                           else
                                for(unsigned int n = 0 ; n < ns_rhs ; ++n){
                                    std::string ptr_name = mat->name() + "_ptr_" + to_string(n);
                                    kss << "__global " << mat->aligned_scalartype() << " * " << ptr_name << " = " << mat->name() << " +  " ;
                                    if(mat->is_transposed())  kss << mat->offset(offset_n + " +  get_group_id(1)*" + to_string(nl_rhs), to_string(n));
                                    else kss << mat->offset("0",offset_n + " +  get_group_id(1)*" + to_string(nl_rhs) + " + " + to_string(n));
                                    kss << ";" << std::endl;
                                    mat->access_name(n,"*" + ptr_name);
                                }
                        }

                    private:
                        utils::kernel_generation_stream & kss;
                        unsigned int ks_rhs;
                        unsigned int ns_rhs;
                        unsigned int nl_rhs;
                        std::string const & offset_n;
                    };


                    struct declare_lhs_global_ptr{

                        declare_lhs_global_ptr(utils::kernel_generation_stream & _kss,
                                               unsigned int _ms_lhs, unsigned int _ks_lhs,
                                               unsigned int _ml_lhs, std::string const & _offset_m) : kss(_kss), ms_lhs(_ms_lhs), ks_lhs(_ks_lhs), ml_lhs(_ml_lhs), offset_m(_offset_m)
                        { }

                        void operator()( mat_infos_base * mat) {
                            if(mat->is_rowmajor()){
                                for(unsigned int m=0; m<ms_lhs; ++m){
                                    std::string ptr_name = mat->name() + "_ptr_" + to_string(m);
                                    kss << "__global " << mat->aligned_scalartype() << " * " << ptr_name << " = " << mat->name() << " + ";
                                    if(mat->is_transposed()) kss << mat->offset(to_string(m),"get_group_id(0)*" + to_string(ml_lhs) + "+" + offset_m );
                                    else kss << mat->offset("get_group_id(0)*" + to_string(ml_lhs) + "+" + offset_m + "+" + to_string(m),"0");
                                    kss << ";" << std::endl;
                                    mat->access_name(m,"*" + ptr_name);
                                }
                            }
                            else{
                                for(unsigned int k=0; k<ks_lhs; ++k){
                                    std::string ptr_name = mat->name() + "_ptr_" + to_string(k);
                                    kss << "__global " << mat->aligned_scalartype() << " * " << ptr_name << " = " << mat->name() << " + " ;
                                    if(mat->is_transposed()) kss << mat->offset("0", to_string(k) + "+" + "get_group_id(0)*" + to_string(ml_lhs) + "+" + offset_m );
                                    else kss << mat->offset( "get_group_id(0)*" + to_string(ml_lhs) + "+" + offset_m, to_string(k));
                                    kss << ";" << std::endl;
                                    mat->access_name(k,"*" + ptr_name);
                                }
                            }
                        }

                    private:
                        utils::kernel_generation_stream & kss;
                        unsigned int ms_lhs;
                        unsigned int ks_lhs;
                        unsigned int ml_lhs;
                        std::string const & offset_m;
                    };

                    struct update_rhs_global_ptr{

                        update_rhs_global_ptr(utils::kernel_generation_stream & _kss, unsigned int _ks, unsigned int _ns_rhs, unsigned int _ks_rhs
                                              ,std::string const & _internal_size1_rhs,
                                              std::string const & _internal_size2_rhs) : kss(_kss), ks(_ks), ns_rhs(_ns_rhs), ks_rhs(_ks_rhs), internal_size1_rhs(_internal_size1_rhs), internal_size2_rhs(_internal_size2_rhs){ }

                        void operator()(mat_infos_base * mat){
                            if(mat->is_rowmajor() && !mat->is_transposed())
                                for(unsigned int k=0 ; k<ks ; ++k)
                                    kss << mat->name() << "_ptr_" << k << " += " << ks_rhs << "*" << internal_size2_rhs << " - " << ns_rhs << ";" << std::endl;
                            else if(mat->is_transposed() && !mat->is_rowmajor())
                                for(unsigned int n=0 ; n<ns_rhs ; ++n)
                                    kss << mat->name() << "_ptr_" << n << " += " << ns_rhs << "*" << internal_size1_rhs << " - " << ks_rhs << ";" << std::endl;
                        }

                    private:
                        utils::kernel_generation_stream & kss;
                        unsigned int ks;
                        unsigned int ns_rhs;
                        unsigned int ks_rhs;
                        std::string const & internal_size1_rhs;
                        std::string const & internal_size2_rhs;
                    };

                    struct update_lhs_global_ptr{

                        update_lhs_global_ptr(utils::kernel_generation_stream & _kss, unsigned int _ks, unsigned int _ms_lhs, unsigned int _ks_lhs
                                              ,std::string const & _internal_size1_lhs,
                                              std::string const & _internal_size2_lhs) : kss(_kss), ks(_ks), ms_lhs(_ms_lhs), ks_lhs(_ks_lhs), internal_size1_lhs(_internal_size1_lhs), internal_size2_lhs(_internal_size2_lhs){ }


                        void operator()(mat_infos_base * mat){
                            if(mat->is_transposed() && mat->is_rowmajor())
                                for(unsigned int m=0 ; m<ms_lhs ; ++m)
                                    kss << mat->name() << "_ptr_" << m << " += " << ks << "*" << internal_size2_lhs << " - " <<  ks_lhs << ";" << std::endl;
                            else if(!mat->is_transposed() && !mat->is_rowmajor())
                                for(unsigned int k=0 ; k<ks_lhs ; ++k)
                                    kss << mat->name() << "_ptr_" << k << " += " << ks_lhs << "*" << internal_size1_lhs << " - " << ms_lhs << ";" << std::endl;
                        }

                    private:
                        utils::kernel_generation_stream & kss;
                        unsigned int ks;
                        unsigned int ms_lhs;
                        unsigned int ks_lhs;
                        std::string const & internal_size1_lhs;
                        std::string const & internal_size2_lhs;
                    };



                    static std::string helper_variable(utils::kernel_generation_stream & kss
                                                        , bool store_in_register
                                                        , std::string const & type
                                                        , std::string const & name
                                                        , std::string const & expr){
                        if(!store_in_register)
                            return expr;
                        kss << type << " " << name << " = " << expr << ";" << std::endl;
                        return name;
                    }

                    template<class MatContainerT>
                    static void fetch_to_local_mem(utils::kernel_generation_stream & kss,
                                                   local_memory<2> const & lmem,
                                                   std::string const & offset,
                                                   unsigned int bound1,
                                                   unsigned int bound2,
                                                   infos_base const & mat_expression,
                                                   MatContainerT & matrices){
                        unsigned int alignment = (*matrices.begin())->alignment();
                        std::string aligned_scalartype = (*matrices.begin())->aligned_scalartype();
                        std::string scalartype = (*matrices.begin())->scalartype();
                        std::string internal_size2 = (*matrices.begin())->internal_size2();
                        std::string internal_size1 = (*matrices.begin())->internal_size1();
                        kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
                        kss << "for(unsigned int i = get_local_id(0)" << " ; i < " << bound1 << "; i+= get_local_size(0)){" << std::endl;
                        kss.inc_tab();
                        kss << "for(unsigned int j = get_local_id(1)" << " ; j < " << bound2 << "; j+= get_local_size(1)){" << std::endl;
                        kss.inc_tab();
                        if((*matrices.begin())->is_rowmajor()){
                            for(typename MatContainerT::iterator it = matrices.begin() ; it!=matrices.end(); ++it){
                                 (*it)->access_name(0,(*it)->name() +  "[" + offset + " + j  + " + internal_size2 + "*i]");
                            }
                            kss << aligned_scalartype << " val = " << mat_expression.generate(0) << ";" << std::endl;
                            kss << "__local " << scalartype << "* ptr = " << lmem.name() << " + i*" << lmem.size2() << "+j*" << alignment<<";" <<std::endl;
                            for(unsigned int a = 0 ; a < alignment ; ++a){
                                if(alignment>1)
                                    kss << "*ptr++ =  val.s" << a << ";" << std::endl;
                                else
                                    kss << "*ptr++ =  val;" << std::endl;
                            }
                        }
                        else{
                            for(typename MatContainerT::iterator it = matrices.begin() ; it!=matrices.end(); ++it){
                                 (*it)->access_name(0,(*it)->name() + "[" + offset + "+ j*" + internal_size1 + " + i]");
                            }
                            kss << aligned_scalartype << " val = " << mat_expression.generate(0) << ";" << std::endl;
                            kss << "__local " << scalartype << "* ptr = " << lmem.name() << " + i*" << alignment * lmem.size2() << "+ j;" <<std::endl;
                            for(unsigned int a = 0 ; a < alignment ; ++a){
                                if(alignment>1)
                                    kss << "*ptr =  val.s" << a << ";" << std::endl;
                                else
                                    kss << "*ptr =  val;" << std::endl;
                                kss << "ptr += " << lmem.size2() << ";" << std::endl;
                            }
                        }

                        kss.dec_tab();
                        kss << "}" << std::endl;
                        kss.dec_tab();
                        kss << "}" << std::endl;
                        kss << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;

                    }

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

                        matmat_prod_infos_base * first_prod = static_cast<matmat_prod_infos_base*>(*matmat_prods_.begin());
                        mat_infos_base* first_lhs = *lhss.begin();
                        mat_infos_base* first_rhs = *rhss.begin();
                        mat_infos_base* first_assigned = assigned.front();

                        bool use_LHS_shared = optimization_profile_->use_LHS_shared();
                        bool use_RHS_shared = optimization_profile_->use_RHS_shared();
                        unsigned int alignment = optimization_profile_->alignment();
                        unsigned int kl = optimization_profile_->kl();
                        unsigned int ks = optimization_profile_->ks();
                        unsigned int ml = optimization_profile_->ml();
                        unsigned int ms = optimization_profile_->ms();
                        unsigned int nl = optimization_profile_->nl();
                        unsigned int ns = optimization_profile_->ns();

                        bool is_lhs_rowmajor = first_lhs->is_rowmajor();
                        bool is_rhs_rowmajor = first_rhs->is_rowmajor();
                        bool is_result_rowmajor = first_assigned->is_rowmajor();

                        bool is_lhs_transposed = first_lhs->is_transposed();
                        bool is_rhs_transposed = first_rhs->is_transposed();

                        std::string lhs_value_scalartype;
                        if(use_LHS_shared) lhs_value_scalartype = first_lhs->scalartype();
                        else lhs_value_scalartype = first_lhs->aligned_scalartype();

                        std::string rhs_value_scalartype;
                        if(use_RHS_shared) rhs_value_scalartype = first_rhs->scalartype();
                        else rhs_value_scalartype = first_rhs->aligned_scalartype();

                        unsigned int ml_res = ml, nl_res = nl, ms_res = ms, ns_res = ns;
                        unsigned int ml_lhs = ml, kl_lhs = kl, ms_lhs = ms, ks_lhs = ks;
                        unsigned int kl_rhs = kl, nl_rhs = nl, ks_rhs = ks, ns_rhs = ns;

                        transform_block(*first_assigned,false,ml_res,nl_res,ms_res,ns_res);
                        transform_block(*first_lhs,use_LHS_shared,ml_lhs,kl_lhs,ms_lhs,ks_lhs);
                        transform_block(*first_rhs,use_RHS_shared,kl_rhs,nl_rhs,ks_rhs,ns_rhs);



                        std::string internal_size1_lhs = first_lhs->internal_size1();
                        std::string internal_size2_lhs = first_lhs->internal_size2();

                        std::string internal_size1_rhs = first_rhs->internal_size1();
                        std::string internal_size2_rhs = first_rhs->internal_size2();

                        std::string internal_size1_res = first_assigned->internal_size1();
                        std::string internal_size2_res = first_assigned->internal_size2();

                        unsigned int lhs_size1 = ml, lhs_size2 = kl;
                        unsigned int rhs_size1 = kl, rhs_size2 = nl;
                        if(is_lhs_transposed) std::swap(lhs_size1, lhs_size2);
                        if(is_rhs_transposed) std::swap(rhs_size1, rhs_size2);

                        local_memory<2> lmem_lhs("local_lhs",lhs_size1,lhs_size2+1,first_lhs->scalartype());
                        local_memory<2> lmem_rhs("local_rhs",rhs_size1,rhs_size2+1,first_lhs->scalartype());

                        //Declaration of results registers
//                        std::string res_table_name(first_prod->repr() + "_res");
                        for(unsigned int m=0; m< ms_res; ++m)
                            for(unsigned int n=0; n < ns_res ; ++n)
                                kss << first_assigned->aligned_scalartype() << " " << first_prod->val_name(m,n) << " = (" << first_assigned->aligned_scalartype() << ")(0) ;" << std::endl;

                        //Declaration of local memories
                        if(use_LHS_shared) kss << lmem_lhs.declare() << ";" << std::endl;
                        if(use_RHS_shared) kss << lmem_rhs.declare() << ";" << std::endl;

                        //Declaration of helpers
                        std::for_each(lhss.begin(),lhss.end(),std::bind1st(std::ptr_fun(transform_size),kss));
                        std::for_each(rhss.begin(),rhss.end(),std::bind1st(std::ptr_fun(transform_size),kss));
                        std::for_each(assigned.begin(),assigned.end(),std::bind1st(std::ptr_fun(transform_size),kss));
                        std::string offset_m = helper_variable(kss,false,"unsigned int", "offset_m", "get_local_id(0)*" + to_string(ms_lhs));
                        std::string offset_n = helper_variable(kss,false,"unsigned int", "offset_n", "get_local_id(1)*" + to_string(ns_rhs));
                        std::string block_num = helper_variable(kss,true,"unsigned int", "block_num", (is_lhs_transposed?internal_size1_lhs:internal_size2_lhs) + '/' + to_string(kl_lhs));

                        //Declaration of pointers and/or offsets to result, rhs, lhs.
                        kss << "__global " << first_assigned->aligned_scalartype() << "* res_ptr = " <<  first_assigned->name() << " + " << first_assigned->offset("get_global_id(0)*" + to_string(ms_res), "get_global_id(1)*" + to_string(ns_res)) << ";" << std::endl;

                        if(use_RHS_shared){
                            if(is_rhs_transposed) kss << "unsigned int offsetRHS = " << first_rhs->offset(" get_group_id(1)*" + to_string(nl_rhs),"0") << ";" << std::endl;
                            else kss << "unsigned int offsetRHS = " << first_rhs->offset("0", " get_group_id(1)*" + to_string(nl_rhs)) << ";" << std::endl;
                        }
                        else{
                            if(is_rhs_transposed)
                                std::for_each(rhss.begin(),rhss.end(),declare_rhs_global_ptr(kss,ns_rhs,ks_rhs,nl_rhs,offset_n));
                            else
                                std::for_each(rhss.begin(),rhss.end(),declare_rhs_global_ptr(kss,ks_rhs,ns_rhs,nl_rhs,offset_n));
                        }

                        if(use_LHS_shared){
                            if(is_lhs_transposed) kss << "unsigned int offsetLHS = " << first_lhs->offset("0", "get_group_id(0)*" + to_string(ml_lhs)) << ";" << std::endl;
                            else kss << "unsigned int offsetLHS = " << first_lhs->offset("get_group_id(0)*" + to_string(ml_lhs), "0") << ";" << std::endl;
                        }
                        else{
                            if(is_lhs_transposed)
                                std::for_each(lhss.begin(),lhss.end(),declare_lhs_global_ptr(kss,ks_lhs,ms_lhs,ml_lhs,offset_m));
                            else
                                std::for_each(lhss.begin(),lhss.end(),declare_lhs_global_ptr(kss,ms_lhs,ks_lhs,ml_lhs,offset_m));
                        }



                        //Main loop
                        kss << "for(unsigned int bl=0 ; bl<" << block_num << " ; ++bl){" << std::endl;
                        kss.inc_tab();

                        //Fetches to local memory if necessary and declares pointers to local memory
                        if(use_LHS_shared){
                            if(is_lhs_transposed) fetch_to_local_mem(kss,lmem_lhs,"offsetLHS",kl_lhs,ml_lhs,first_prod->lhs(),lhss);
                            else fetch_to_local_mem(kss,lmem_lhs,"offsetLHS",ml_lhs,kl_lhs,first_prod->lhs(),lhss);
                            unsigned int upper_bound = is_lhs_transposed?ks_lhs:ms_lhs;
                            for(unsigned int m=0; m<upper_bound; ++m){
                                 kss << "__local " << lhs_value_scalartype << "* ptr_lhs_" << m << " = local_lhs + " ;
                                if(is_lhs_transposed) kss << m*lmem_lhs.size2() << " + " << offset_m ;
                                else kss << "(" << offset_m << "+" << m << ")" << "*" << lmem_lhs.size2() ;
                                kss << ";" << std::endl;
                            }
                        }

                        if(use_RHS_shared){
                            if(is_rhs_transposed) fetch_to_local_mem(kss,lmem_rhs,"offsetRHS",nl_rhs,kl_rhs,first_prod->rhs(),rhss);
                            else fetch_to_local_mem(kss,lmem_rhs,"offsetRHS",kl_rhs,nl_rhs,first_prod->rhs(),rhss);
                            unsigned int upper_bound = is_rhs_transposed?ns_rhs:ks_rhs;
                            for(unsigned int k=0; k<upper_bound; ++k){
                                kss << "__local " << rhs_value_scalartype << "* ptr_rhs_" << k << " = local_rhs + " ;
                                if(is_rhs_transposed) kss << "(" << offset_n << "+" << k << ")*" << lmem_rhs.size2();
                                else kss << k*lmem_rhs.size2() << " + " << offset_n;
                                kss << ";" << std::endl;
                            }
                        }


                        kss << " for(unsigned int bs=0 ; bs < " << kl/ks  << " ; ++bs){" << std::endl;
                        kss.inc_tab();


                        unsigned int upperbound_1_rhs = is_rhs_transposed?ns_rhs:ks_rhs;
                        unsigned int upperbound_2_rhs = is_rhs_transposed?ks_rhs:ns_rhs;
                        for(unsigned int k = 0 ; k < upperbound_1_rhs ; ++k){
                            for(unsigned int n=0 ; n < upperbound_2_rhs ; ++n){
                                kss << rhs_value_scalartype << " val_rhs_" << k << "_" << n << " = " ;
                                if(use_RHS_shared ) kss << "* ptr_rhs_" << k << "++";
                                else{
                                    if(is_rhs_rowmajor) kss << first_prod->rhs().generate(k);
                                    else kss  << first_prod->rhs().generate(n);
                                }
                                kss << ";";
                                if( !use_RHS_shared ){
                                    for(std::set<mat_infos_base*,viennacl::generator::deref_less>::iterator it = rhss.begin() ; it!=rhss.end() ; ++it){
                                        if(is_rhs_rowmajor)kss << "++" << (*it)->name() << "_ptr_" << k << ";" ;
                                        else kss << "++" << (*it)->name() << "_ptr_" << n << ";" ;
                                    }
                                }
                                kss << std::endl;
                            }
                        }



                        unsigned int upperbound_1_lhs = is_lhs_transposed?ms_lhs:ks_lhs;
                        unsigned int upperbound_2_lhs = is_lhs_transposed?ks_lhs:ms_lhs;
                        for(unsigned int k = 0 ; k < upperbound_1_lhs ; ++k){
                            for(unsigned int m=0 ; m < upperbound_2_lhs ; ++m){
                                kss << lhs_value_scalartype << " " << "val_lhs_" << m << "_" << k << " = ";
                                if(use_LHS_shared) kss <<  "* ptr_lhs_" << m << "++" ;
                                else if(is_lhs_rowmajor) kss << first_prod->lhs().generate(m);
                                else kss << first_prod->lhs().generate(k);
                                kss << ";";
                                if( !use_LHS_shared ){
                                    for(std::set<mat_infos_base*,viennacl::generator::deref_less>::iterator it = lhss.begin() ; it!=lhss.end() ; ++it){
                                        if(is_lhs_rowmajor) kss << "++" << (*it)->name() << "_ptr_" << m << ";" ;
                                        else kss << "++" << (*it)->name() << "_ptr_" << k << ";" ;
                                    }
                                }
                                kss << std::endl;
                            }
                        }



                        for(unsigned int k = 0 ; k < ks ; ++k){
                            for(unsigned int n=0 ; n < ns_res ; ++n){
                                for(unsigned int m=0 ; m < ms_res ; ++m){
                                    for(unsigned int a = 0; a<alignment; ++a){

                                        int ind_lhs_1 = m;
                                        int ind_lhs_2 = k;
                                        int ind_s_lhs=a;

                                        int ind_rhs_1=k;
                                        int ind_rhs_2=n;
                                        int ind_s_rhs=a;

                                        bool is_vectorized_lhs = false;
                                        bool is_vectorized_rhs = false;

                                        if(is_result_rowmajor){
                                            if(is_lhs_transposed) std::swap(ind_lhs_1,ind_lhs_2);

                                            if(!use_LHS_shared){
                                                if(is_lhs_rowmajor){
                                                    ind_s_lhs = ind_lhs_2%alignment;
                                                    ind_lhs_2 /= alignment;
                                                }
                                                else{
                                                    ind_s_lhs = ind_lhs_1%alignment;
                                                    ind_lhs_1 /= alignment;
                                                }
                                            }
                                        }
                                        else{
                                            if(use_LHS_shared){
                                                ind_lhs_1 = ind_lhs_1*alignment+a;
                                            }
                                            else{
                                                if((is_lhs_rowmajor && !is_lhs_transposed)
                                                        ||(!is_lhs_rowmajor && is_lhs_transposed)){
                                                    ind_lhs_1 = ind_lhs_1*alignment+a;
                                                    ind_s_lhs = ind_lhs_2%alignment;
                                                    ind_lhs_2 /= alignment;

                                                }
                                            }
                                            if(is_lhs_transposed) std::swap(ind_lhs_1,ind_lhs_2);
                                        }

                                        if(is_result_rowmajor){
                                            if(use_RHS_shared){
                                                ind_rhs_2 = ind_rhs_2*alignment+a;
                                            }
                                            else{
                                                if((!is_rhs_rowmajor && !is_rhs_transposed)
                                                    ||(is_rhs_rowmajor && is_rhs_transposed)){
                                                    ind_rhs_2 = ind_rhs_2*alignment+a;
                                                    ind_s_rhs = ind_rhs_1%alignment;
                                                    ind_rhs_1 = ind_rhs_1/alignment;
                                                }
                                                else if( (is_rhs_rowmajor && !is_rhs_transposed) ){
                                                    is_vectorized_rhs=true;
                                                }
                                            }
                                            if(is_rhs_transposed) std::swap(ind_rhs_1,ind_rhs_2);
                                        }
                                        else{
                                            if(is_rhs_transposed) std::swap(ind_rhs_1,ind_rhs_2);
                                            if(!use_RHS_shared){
                                                if(is_rhs_rowmajor){
                                                    ind_s_rhs = ind_rhs_2%alignment;
                                                    ind_rhs_2/=alignment;
                                                }
                                                else{
                                                    ind_s_rhs = ind_rhs_1%alignment;
                                                    ind_rhs_1/=alignment;
                                                }
                                            }
                                        }

                                        bool is_vectorized = is_vectorized_lhs || is_vectorized_rhs;

                                        std::ostringstream res_oss;
                                        std::ostringstream lhs_oss;
                                        std::ostringstream rhs_oss;

                                        res_oss << first_prod->val_name(m,n);
                                        if(!is_vectorized && alignment>1) res_oss << ".s" << a;

                                        lhs_oss << "val_lhs_" << ind_lhs_1 << "_" << ind_lhs_2;
                                        if(!is_vectorized_lhs && !use_LHS_shared && alignment>1) lhs_oss << ".s" << ind_s_lhs;


                                        rhs_oss << "val_rhs_" << ind_rhs_1 << "_" << ind_rhs_2;
                                        if(!is_vectorized_rhs && !use_RHS_shared && alignment>1) rhs_oss << ".s" << ind_s_rhs;

                                        kss << first_prod->update_val(res_oss.str(),lhs_oss.str(), rhs_oss.str()) << ";" << std::endl;


                                        if(is_vectorized)
                                            break;
                                    }
                                }
                            }
                        }


                        if(use_RHS_shared){
                            for(unsigned int k=0 ; k<ks ; ++k)
                                if(!is_rhs_transposed) kss << "ptr_rhs_" << k << " += " << ks_rhs*lmem_rhs.size2() - ns_rhs << ";" << std::endl;
                        }
                        else{
                            if(is_rhs_transposed)
                                std::for_each(rhss.begin(),rhss.end(),update_rhs_global_ptr(kss,ks,ks_rhs,ns_rhs,internal_size1_rhs,internal_size2_rhs));
                            else
                                 std::for_each(rhss.begin(),rhss.end(),update_rhs_global_ptr(kss,ks,ns_rhs,ks_rhs,internal_size1_rhs,internal_size2_rhs));
                        }



                        if(use_LHS_shared){
                            for(unsigned int m=0 ; m<ks_lhs ; ++m)
                                if(is_lhs_transposed) kss << "ptr_lhs_" << m << " += " << ks*lmem_lhs.size2() - ms_lhs << ";" << std::endl;
                        }
                        else{
                            if(is_lhs_transposed)
                                std::for_each(lhss.begin(),lhss.end(),update_lhs_global_ptr(kss,ks,ks_lhs,ms_lhs,internal_size1_lhs,internal_size2_lhs));
                            else
                                std::for_each(lhss.begin(),lhss.end(),update_lhs_global_ptr(kss,ks,ms_lhs,ks_lhs,internal_size1_lhs,internal_size2_lhs));
                        }



                        kss.dec_tab();
                        kss << "}" << std::endl;

                        if(use_LHS_shared){
                            if(is_lhs_transposed){
                                if(is_lhs_rowmajor)
                                    kss << "offsetLHS += " << kl_lhs << "*" << internal_size2_lhs << ";" << std::endl;
                                else
                                    kss << "offsetLHS += " << kl_lhs  << ";" << std::endl;
                            }
                            else{
                                if(is_lhs_rowmajor)
                                    kss << "offsetLHS += " << kl_lhs << ";" << std::endl;
                                else
                                    kss << "offsetLHS += " << kl_lhs << "*" << internal_size1_lhs << ";" << std::endl;
                            }

                        }

                        if(use_RHS_shared){
                            if(is_rhs_transposed){
                                if(is_rhs_rowmajor)
                                    kss << "offsetRHS += " << kl_rhs << ";" << std::endl;
                                else
                                    kss << "offsetRHS += " << kl_rhs << "*" << internal_size1_rhs << ";" << std::endl;
                            }
                            else{
                                if(is_rhs_rowmajor)
                                    kss << "offsetRHS += " << kl_rhs << "*" << internal_size2_rhs << ";" << std::endl;
                                else
                                    kss << "offsetRHS += " << kl_rhs << ";" << std::endl;
                            }
                        }

                        kss.dec_tab();
                        kss << "}" << std::endl;

                        if(first_assigned->is_rowmajor()){
                            for(unsigned int m=0 ; m < ms_res ; ++m){
                                for(unsigned int n=0 ; n < ns_res ; ++n){
                                    kss << "*res_ptr++=" << first_prod->val_name(m,n) << ";" << std::endl;
                                }
                                kss << "res_ptr+=" << internal_size2_res << " - " << ns_res << ";" << std::endl;
                            }
                        }
                        else{
                            for(unsigned int n=0 ; n < ns_res ; ++n){
                                for(unsigned int m=0 ; m < ms_res ; ++m){
                                    kss << "*res_ptr++=" << first_prod->val_name(m,n) << ";" << std::endl;
                                }
                                kss << "res_ptr+=" << internal_size1_res << " - " << ms_res << ";" << std::endl;
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
