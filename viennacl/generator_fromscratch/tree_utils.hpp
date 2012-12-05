#ifndef VIENNACL_GENERATOR_TREE_UTILS_HPP
#define VIENNACL_GENERATOR_TREE_UTILS_HPP

#include "list"
#include <algorithm>
#include "viennacl/generator_fromscratch/symbolic_types.hpp"

namespace viennacl{

    namespace generator{

        struct tree_utils{

            typedef std::list<infos_base*>  tree_list_t;
            typedef std::list<infos_base*> leaves_list_t;

            template<class Pred>
            static leaves_list_t extract_if(infos_base* const node, Pred pred, bool inspect_binary_leaves=true){
                leaves_list_t res;
                if(arithmetic_tree_infos_base * p = dynamic_cast<arithmetic_tree_infos_base *>(node)){
                        leaves_list_t tmplhs(extract_if(&p->lhs(),pred));
                        leaves_list_t tmprhs(extract_if(&p->rhs(),pred));
                        res.merge(tmplhs);
                        res.merge(tmprhs);
                }
                else if(unary_tree_infos_base * p = dynamic_cast<unary_tree_infos_base *>(node)){
                    leaves_list_t tmp(extract_if(&p->sub(),pred));
                    res.merge(tmp);
                }
                if(pred(node)){
                    res.push_back(node);
                }
                return res;
            }

            template<class Pred>
            static leaves_list_t extract_if(tree_list_t const & trees, Pred pred, bool inspect_binary_leaves=true){
                leaves_list_t res;
                for(tree_list_t::const_iterator it = trees.begin() ; it != trees.end() ; ++it){
                    leaves_list_t tmp(extract_if(*it,pred));
                    res.merge(tmp);
                }
                return res;
            }

            template<class T>
            static std::list<T *> extract_type(tree_list_t const & trees, bool inspect_binary_leaves=true){
                leaves_list_t tmp(tree_utils::extract_if(trees,is_type<kernel_argument>));
                std::list<T*> res(tmp.size());
                std::transform(tmp.begin(), tmp.end(), res.begin(),Base2Target<infos_base,T>());
                return res;
            }
        };
    }
}
#endif // GENERATE_HEADER_HPP
