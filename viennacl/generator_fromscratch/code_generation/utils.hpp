#ifndef VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP

#include <list>
#include <set>

namespace viennacl{

    namespace generator{


        namespace code_generation{


            namespace utils{


                template<class T>
                static bool is_type(infos_base* p){
                    return dynamic_cast<T *>(p);
                }

                template<class T>
                static bool is_pointed_value_eq(T* a, T* b){
                    return (*a == *b);
                }

                template<class T>
                static bool is_pointed_value_inf(T* a, T* b){
                    return (*a) <= (*b);
                }

                template<class Pred>
                static std::list<infos_base*> extract_if(infos_base* const node, Pred pred, bool inspect_nested_leaves=true){
                    std::list<infos_base*> res;

                    if(inspect_nested_leaves){
                        if(arithmetic_tree_infos_base * p = dynamic_cast<arithmetic_tree_infos_base *>(node)){
                                std::list<infos_base*> tmplhs(extract_if(&p->lhs(),pred));
                                std::list<infos_base*> tmprhs(extract_if(&p->rhs(),pred));
                                res.merge(tmplhs);
                                res.merge(tmprhs);
                        }
                        else if(unary_tree_infos_base * p = dynamic_cast<unary_tree_infos_base *>(node)){
                            std::list<infos_base*> tmp(extract_if(&p->sub(),pred));
                            res.merge(tmp);
                        }
                    }
                    if(pred(node)){
                        res.push_back(node);
                    }
                    return res;
                }

                template<class Pred>
                static std::list<infos_base*> extract_if(std::list<infos_base*> const & trees, Pred pred, bool inspect_nested_leaves=true){
                    std::list<infos_base*> res;
                    for(std::list<infos_base*>::const_iterator it = trees.begin() ; it != trees.end() ; ++it){
                        std::list<infos_base*> tmp(extract_if(*it,pred,inspect_nested_leaves));
                        res.merge(tmp);
                    }
                    return res;
                }


                template<class T, class B>
                static std::list<T *> cast(std::list<B *> const & in){
                    std::list<T*> res(in.size());
                    std::transform(in.begin(), in.end(), res.begin(),Base2Target<B,T>());
                    return res;
                }

                template<class T>
                static std::list<T *> extract_cast(std::list<infos_base*> const & trees){
                    return cast<T>(extract_if(trees,is_type<T>));
                }

                template<class T>
                class cache_manager{
                public:
                    cache_manager(std::list<T * > const & expressions,  std::ostringstream & oss, std::set<T *>& cached_entries) : expressions_(expressions)
                                                                                                                                               , oss_(oss)
                      ,cached_entries_(cached_entries){         }

                    void fetch_entries(std::string const & idx){
                        for(typename std::list<T * >::iterator it = expressions_.begin() ; it != expressions_.end() ; ++it){
                            T * p = *it;
                            if(cached_entries_.insert(p).second){
                                p->access_name(p->name()+"_val");
                                oss_ << p->scalartype() << " " << p->generate() << " = " << p->name() << "[" << idx << "];\n";
                            }
                        }
                    }

                    void writeback_entries(std::string const & idx){
                        for(typename std::list<T * >::iterator it = expressions_.begin() ; it != expressions_.end() ; ++it){
                            T * p = *it;
                            if(p->is_modified())
                                oss_<< p->name() << "[" << idx << "]"<< " = "  << p->generate() << ";\n";
                        }
                    }

                private:
                    std::list<T * > expressions_;
                    std::ostringstream & oss_;
                    std::set<T *> & cached_entries_;

                };

            }

        }

    }

}
#endif // UTILS_HPP
