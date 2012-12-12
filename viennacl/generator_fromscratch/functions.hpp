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

    class symbolic_function : public function_base{
    public:
        template<class T1>
        symbolic_function(std::string const & name
                        ,std::string const & expr
                          ,T1 const & t1) : function_base(name), expr_(expr){
            args_map_.insert(std::make_pair("_1_", new T1(t1)));
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
        unsigned int n_args_;
    };





}

}
#endif // FUNCTIONS_HPP
