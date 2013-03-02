#ifndef VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_FRONTEND_HPP

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <typeinfo>

#include "viennacl/generator/symbolic_types_base.hpp"
#include "viennacl/generator/code_generation/backend.hpp"
#include "viennacl/generator/code_generation/utils.hpp"

#include "viennacl/tools/shared_ptr.hpp"

#ifdef VIENNACL_ENABLE_AUTOTUNE
#include "viennacl/io/kernel_parameters.hpp"
#endif

namespace viennacl{

    namespace generator{


        template<class ArgumentsT>
        void set_arguments(viennacl::ocl::kernel & k, ArgumentsT const & args){
            unsigned int counter=0;
            for(typename ArgumentsT::const_iterator iit = args.begin(); iit != args.end() ; ++iit){
                (*iit)->enqueue(counter,k);
            }
        }




        namespace code_generation{

            enum kernel_type_t{
                BLAS1_TYPE,
                BLAS2_TYPE,
                BLAS3_TYPE
            };

            template<class T>
            struct kernel_type_of;

            template<>
            struct kernel_type_of<blas3_optimization_profile>{
                static const kernel_type_t value = BLAS3_TYPE;
            };

            template<>
            struct kernel_type_of<blas1_optimization_profile>{
                static const kernel_type_t value = BLAS1_TYPE;
            };

            class kernel_infos_t{
            public:

                template<class ProfType>
                kernel_infos_t(infos_base *op, ProfType const & prof_model) : type_(kernel_type_of<ProfType>::value), optimization_profile_(new ProfType(prof_model)){
                    trees_.push_back(op);
                }


                typedef std::list<kernel_argument*> arguments_t;

                kernel_type_t type(){ return type_; }

                std::list<infos_base*> & trees(){ return trees_; }

                arguments_t & arguments(){ return arguments_; }

                code_generation::optimization_profile* profile() { return optimization_profile_.get(); }

                void config_nd_range(viennacl::ocl::kernel & k) const{
                    if(type_==BLAS3_TYPE){
                        matrix_expression_infos_base * p = static_cast<matrix_expression_infos_base*>(trees_.front());
                        mat_infos_base * mat = static_cast<mat_infos_base*>(&p->lhs());
                        blas3_optimization_profile* prof = static_cast<blas3_optimization_profile*>(optimization_profile_.get());
                        prof->set_global_sizes(mat);
                    }
                    optimization_profile_->config_nd_range(k);
                }

            private:
                arguments_t arguments_;
                std::list<infos_base*> trees_;
                kernel_type_t type_;
                viennacl::tools::shared_ptr<code_generation::optimization_profile> optimization_profile_;
            };



            template<class T>
            static bool find_pred(T* t1, T* t2){
                return t1->handle()==t2->handle();
            }

            template<class T, class Pred>
            static void extract_to_list(infos_base* root, std::list<T*> & args, Pred pred){
                if(arithmetic_tree_infos_base* p = dynamic_cast<arithmetic_tree_infos_base*>(root)){
                    extract_to_list(&p->lhs(), args,pred);
                    extract_to_list(&p->rhs(),args,pred);
                }
                else if(function_base* p = dynamic_cast<function_base*>(root)){
                    std::list<infos_base*> func_args(p->args());
                    for(std::list<infos_base*>::const_iterator it = func_args.begin(); it!= func_args.end(); ++it){
                        extract_to_list(*it,args,pred);
                    }
                }
                else if(inprod_infos_base* p = dynamic_cast<inprod_infos_base*>(root)){
                    if(p->step() == inprod_infos_base::compute){
                        extract_to_list(&p->lhs(), args,pred);
                        extract_to_list(&p->rhs(),args,pred);
                    }
                }
                if(T* t = dynamic_cast<T*>(root)){
                    if(pred(t)){
                        typename std::list<T*>::iterator it = std::find_if(args.begin(),args.end(),std::bind2nd(std::ptr_fun(find_pred<T>),t));
                        if(it==args.end()) args.push_back(t);
                    }
                }
            }


            class kernel_generator{
            private:

                void generate_headers(){
                    kss_ << "__kernel void " + kernel_name_ + "(";
                    for(std::list<infos_base*>::iterator it = kernel_infos_.trees().begin() ; it!= kernel_infos_.trees().end() ; ++it){
                        extract_to_list(*it,kernel_infos_.arguments(),utils::is_type<kernel_argument>());
                    }
                    for(std::list<kernel_argument*>::iterator it=kernel_infos_.arguments().begin(); it!=kernel_infos_.arguments().end();++it){
                        if(it!=kernel_infos_.arguments().begin()) kss_ << ',';
                        kss_ << (*it)->arguments_string() << std::endl ;
                    }
                    kss_ << ")" << std::endl;
                }

                void generate_sources(){
                    kss_<<"{"<< std::endl;
                    kss_.inc_tab();
                    std::list<infos_base *> vec_exprs;
                    std::list<infos_base *> scal_exprs;
                    std::list<infos_base *> mat_exprs;
                    for(std::list<infos_base*>::const_iterator it = kernel_infos_.trees().begin(); it!=kernel_infos_.trees().end();++it){
                        if(utils::is_type<vector_expression_infos_base>()(*it))
                            vec_exprs.push_back(*it);
                        else if(utils::is_type<scalar_expression_infos_base>()(*it)
                                ||utils::is_type<inprod_infos_base>()(*it))
                            scal_exprs.push_back(*it);
                        else
                            mat_exprs.push_back(*it);
                    }
                    if(kernel_infos_.type()==BLAS1_TYPE){
                        code_generation::blas1_generator gen(vec_exprs,mat_exprs,scal_exprs,static_cast<blas1_optimization_profile *>(kernel_infos_.profile()));
                        gen(kss_);
                    }
                    else if(kernel_infos_.type()==BLAS3_TYPE){
                        code_generation::blas3_generator gen(mat_exprs,static_cast<blas3_optimization_profile*>(kernel_infos_.profile()));
                        gen(kss_);
                    }
                    kss_.dec_tab();
                    kss_<<"}"<< std::endl;
                }

            public:
                kernel_generator(kernel_infos_t & kernel_infos
                                 , std::string const & kernel_name
                                 , code_generation::utils::kernel_generation_stream & kss) : kernel_infos_(kernel_infos)
                                                                     , kernel_name_(kernel_name)
                                                                     , kss_(kss){
                    kernel_infos_.profile()->apply(kernel_infos_.trees());
                }

                void generate(){
                    generate_headers();
                    generate_sources();
                }


            private:
                kernel_infos_t & kernel_infos_;
                std::string kernel_name_;
                utils::kernel_generation_stream & kss_;
            };



            class operations_manager{
            private:
                typedef std::list<viennacl::tools::shared_ptr<infos_base> > operations_t;

#ifdef VIENNACL_ENABLE_AUTOTUNE
                template<class T>
                T get_param(viennacl::io::parameter_database & paras, std::string const & root, std::string const & name) const{
                    T res;
                    std::string par_str = root+"/params/param[name='" + name + "']/value/text()";
                    pugi::xpath_node_set par_res = paras.doc.select_nodes(par_str.c_str());
                    std::stringstream ss;
                    ss << par_res[0].node().value();
                    ss >> res;
                    return res;
                }
#endif

                blas3_optimization_profile load_blas3_profile(matrix_expression_infos_base const * operation) const{
                    cl_device_id id = viennacl::ocl::current_device().id();
                    blas3_optimization_profile default_profile = builtin_dabase[std::make_pair(ocl::info<CL_DEVICE_VENDOR_ID>(id), ocl::info<CL_DEVICE_TYPE>(id))][operation->simplified_repr()];

#ifndef VIENNACL_ENABLE_AUTOTUNE
                    return default_profile;
#else

                    const char * tmp = std::getenv("VIENNACL_PARAMS_PATH");
                    if(tmp==NULL){
//                        std::cerr << "Tuner: Environment variable VIENNACL_PARAMS_PATH not set ... Falling back on default" << std::endl;
//                        return default_profile;
                        tmp = "./";
                    }
                    std::string VCL_PARAMS_PATH(tmp);

                    std::string filename(VCL_PARAMS_PATH + "/blas3_parameters.xml");
                    viennacl::io::parameter_database  paras;
                    paras.load(filename);

                    if(paras.doc.empty()){
                        std::cerr << "Tuner : XML Results are empty or nonexistant ... Falling back on default" << std::endl;
                        return default_profile;
                    }

                    std::string devname = viennacl::ocl::current_device().name();
                    // check if tune parameters for the current device are present
                    std::string          device_str = "/parameters/devices/device[contains(name,'"+devname+"')]";
                    pugi::xpath_node_set device_res = paras.doc.select_nodes(device_str.c_str());

                    if(device_res.empty()){
                        std::cerr << "Tuner: There are no parameters existing for the device : " << devname << " ! ... Falling back on default" << std::endl;
                        return default_profile;
                    }
                    else if(device_res.size()>1){
                        std::cerr << "Tuner: Existing multiple definitions for the device : " << devname << " ! ... Falling back on default" << std::endl;
                        return default_profile;
                    }

                    // check if tune parameters for float exist
                    std::string          name_str = device_str+"/kernels/kernel[name='"+operation->simplified_repr()+"']";
                    pugi::xpath_node_set name_res = paras.doc.select_nodes(name_str.c_str());

                    if(name_res.size() == 0){
                        std::cerr << "Tuner: There are no parameters for the kernel : " << operation->simplified_repr() << "! ... Falling back on default" << std::endl;
                        return default_profile;
                    }

                    blas3_optimization_profile res(get_param<unsigned int>(paras,name_str,"ml"),get_param<unsigned int>(paras,name_str,"kl"),get_param<unsigned int>(paras,name_str,"nl"),
                                                      get_param<unsigned int>(paras,name_str,"ms"),get_param<unsigned int>(paras,name_str,"ks"),get_param<unsigned int>(paras,name_str,"ns"),
                                                      get_param<bool>(paras,name_str,"lhs_storage"),get_param<bool>(paras,name_str,"rhs_storage"),get_param<unsigned int>(paras,name_str,"alignment")
                                );

//                    std::cout << operation->repr() << " " << res << std::endl;
                    return res;
#endif
                }

                void init(){
                    if(!kernels_list_.empty()) return;
                    for(typename operations_t::const_iterator it = operations_.begin() ; it!=operations_.end() ; ++it){
                        infos_base* p = it->get();
                        std::list<infos_base*> inprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<inprod_infos_base>()));
                        std::list<infos_base*> matmatprods(utils::filter<utils::EXTRACT_IF>(p,utils::is_type<matmat_prod_infos_base>()));
                        kernel_infos_t* last = NULL;

                        blas1_optimization_profile blas1_model;
                        if(overriden_blas1_model_.get()) blas1_model = *overriden_blas1_model_;

                        if(inprods.size()){
                            assert(matmatprods.size()==0 && "INVALID KERNEL !");
                            if(kernels_list_.size()){
                                if(kernels_list_.back().type() != BLAS1_TYPE)
                                    kernels_list_.push_back(kernel_infos_t(p,blas1_model));
                                else
                                    kernels_list_.back().trees().merge(inprods);
                            }
                            else{
                                kernels_list_.push_back(kernel_infos_t(p,blas1_model));
                            }
                        }
                        else if(matmatprods.size()){
                            blas3_optimization_profile blas3_model = load_blas3_profile(static_cast<matrix_expression_infos_base*>(p));
                            if(overriden_blas3_model_.get()) blas3_model = *overriden_blas3_model_;
                            blas3_optimization_profile prof(blas3_model);
                            if(kernels_list_.size()){
                                if(kernels_list_.back().type() != BLAS3_TYPE)
                                    kernels_list_.push_back(kernel_infos_t(p,prof));
                                else
                                    kernels_list_.back().trees().push_back(p);
                            }
                            else{
                                kernels_list_.push_back(kernel_infos_t(p,prof));
                            }
                        }
                        else{
                            if(kernels_list_.size()){
                                if(kernels_list_.back().type() != BLAS1_TYPE)
                                    kernels_list_.push_back(kernel_infos_t(p,blas1_model));
                                else
                                    kernels_list_.back().trees().push_back(p);
                            }
                            else{
                                kernels_list_.push_back(kernel_infos_t(p,blas1_model));
                            }
                        }

                    }
                }

            public:

                operations_manager(){

                }

                void override_blas1_model(blas1_optimization_profile const & o){
                    overriden_blas1_model_.reset(new blas1_optimization_profile(o));
                }

                void override_blas3_model(blas3_optimization_profile const & o){
                    overriden_blas3_model_.reset(new blas3_optimization_profile(o));
                }

                template<class T>
                void add(T const & op){
                    operations_.push_back(viennacl::tools::shared_ptr<infos_base>(new T(op)));
                }

                void flush(){
                    operations_.clear();
                }

                std::list<kernel_infos_t> get_kernels_list(){
                    init();
                    return kernels_list_;
                }

                std::string repr(){
                    init();
                    std::string res;
                    for(std::list<kernel_infos_t>::iterator it = kernels_list_.begin() ; it !=kernels_list_.end() ; ++it){
                        for(std::list<infos_base*>::iterator iit = it->trees().begin() ; iit != it->trees().end() ; ++iit){
                            res += (*iit)->repr();
                        }
                        res+=it->profile()->repr();
                    }
                    return res;
                }

                std::string get_source_code( std::map<std::string, generator::code_generation::kernel_infos_t> & kernels_infos){
                    init();

                    std::ostringstream oss;
                    code_generation::utils::kernel_generation_stream kss(oss);
                    kss << "#if defined(cl_khr_fp64)\n";
                    kss <<  "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n";
                    kss <<  "#elif defined(cl_amd_fp64)\n";
                    kss <<  "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n";
                    kss <<  "#endif\n";

                    for(std::list<kernel_infos_t>::iterator it = kernels_list_.begin() ; it !=kernels_list_.end() ; ++it){
                        std::string name("_k"+to_string(std::distance(kernels_list_.begin(),it)));
                        kernel_infos_t & infos = kernels_infos.insert(std::make_pair(name,*it)).first->second;
                        kss <<  "__attribute__((reqd_work_group_size(" << infos.profile()->local_work_size(0)
                                                                        << "," << std::max(infos.profile()->local_work_size(1),(unsigned int)1)
                                                                        << ",1)))" << std::endl;
                        code_generation::kernel_generator kg(infos,name,kss);
                        kg.generate() ;
                    }
                    return oss.str();
                }

            private:
                operations_t operations_;
                viennacl::tools::shared_ptr<blas1_optimization_profile> overriden_blas1_model_;
                viennacl::tools::shared_ptr<blas3_optimization_profile> overriden_blas3_model_;
                std::list<kernel_infos_t> kernels_list_;

            };

        }


    }

}
#endif // KERNEL_GENERATOR_FRONTEND_HPP
