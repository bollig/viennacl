#ifndef VIENNACL_GENERATOR_FUNCTIONS_HPP
#define VIENNACL_GENERATOR_FUNCTIONS_HPP


#include "viennacl/generator_fromscratch/symbolic_types_base.hpp"
#include "viennacl/generator_fromscratch/dummy_types.hpp"
#include <map>
#include <list>
#include <cassert>

namespace viennacl{

namespace generator{

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

    class inline_function : public function_base{
    public:
        inline_function(std::string const & name
                        ,std::string const & expr) : function_base(name),expr_(expr){
            n_args_ = 0;
            bool keep_going = true;
            while(keep_going){
                std::string current_arg = "_"+to_string(n_args_+1)+"_";
                if(expr_.find(current_arg)!=std::string::npos)
                    ++n_args_;
                else
                    keep_going=false;
            }
            assert(n_args_>0 && "\nNo argument specified for the function\n"
                                "\nRecall : 1st arg : _1_\n"
                                "\n         2nd arg : _2_\n"
                                          "...");
        }

        std::list<infos_base*> args() const{
            std::list<infos_base*> res;
            for(args_map_t::const_iterator it = args_map_.begin() ; it!= args_map_.end() ; ++it)
                res.push_back(it->second.get());
            return res;
        }

        template<class T1>
        inline_function& operator()(T1 const & t1){
            assert(n_args_==1);
            args_map_.insert(std::make_pair("_1_",viennacl::tools::shared_ptr<infos_base>(new T1(t1))));
            return *this;
        }

        template<class T1, class T2>
        inline_function& operator()(T1 const & t1, T2 const & t2){
            assert(n_args_==2);
            args_map_.insert(std::make_pair("_1_",viennacl::tools::shared_ptr<infos_base>(new T1(t1))));
            args_map_.insert(std::make_pair("_2_",viennacl::tools::shared_ptr<infos_base>(new T2(t2))));
            return *this;
        }

        template<class T1, class T2, class T3>
        inline_function& operator()(T1 const & t1, T2 const & t2, T3 const & t3){
            assert(n_args_==3);
            args_map_.insert(std::make_pair("_1_",viennacl::tools::shared_ptr<infos_base>(new T1(t1))));
            args_map_.insert(std::make_pair("_2_",viennacl::tools::shared_ptr<infos_base>(new T2(t2))));
            args_map_.insert(std::make_pair("_3_",viennacl::tools::shared_ptr<infos_base>(new T3(t3))));
            return *this;
        }


        template<class T1, class T2, class T3, class T4>
        inline_function& operator()(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4){
            assert(n_args_==4);
            args_map_.insert(std::make_pair("_1_",viennacl::tools::shared_ptr<infos_base>(new T1(t1))));
            args_map_.insert(std::make_pair("_2_",viennacl::tools::shared_ptr<infos_base>(new T2(t2))));
            args_map_.insert(std::make_pair("_3_",viennacl::tools::shared_ptr<infos_base>(new T3(t3))));
            args_map_.insert(std::make_pair("_4_",viennacl::tools::shared_ptr<infos_base>(new T4(t4))));
            return *this;
        }

        template<class T1, class T2, class T3, class T4, class T5>
        inline_function& operator()(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5){
            assert(n_args_==4);
            args_map_.insert(std::make_pair("_1_",viennacl::tools::shared_ptr<infos_base>(new T1(t1))));
            args_map_.insert(std::make_pair("_2_",viennacl::tools::shared_ptr<infos_base>(new T2(t2))));
            args_map_.insert(std::make_pair("_3_",viennacl::tools::shared_ptr<infos_base>(new T3(t3))));
            args_map_.insert(std::make_pair("_4_",viennacl::tools::shared_ptr<infos_base>(new T4(t4))));
            args_map_.insert(std::make_pair("_4_",viennacl::tools::shared_ptr<infos_base>(new T5(t5))));
            return *this;
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
        unsigned int n_args_;
    };





}

}
#endif // FUNCTIONS_HPP
