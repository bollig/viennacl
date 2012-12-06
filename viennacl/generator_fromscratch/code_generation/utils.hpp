#ifndef VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP

#include <list>
#include <set>
#include <algorithm>

namespace viennacl{

    namespace generator{


        namespace code_generation{


            namespace utils{


                template<class T>
                static bool is_type(infos_base* p){
                    return dynamic_cast<T *>(p);
                }

                struct less {
                  template<class T>
                  bool operator()(T &a, T &b) {
                    return std::less<T>()(a, b);
                  }
                };

                struct deref_less {
                  template<class T>
                  bool operator()(T a, T b) {
                    return less()(*a, *b);
                  }
                };

                struct double_deref_less {
                  template<class T>
                  bool operator()(T a, T b) {
                    return less()(**a, **b);
                  }
                };

                template<class T>
                struct deref_t{ typedef deref_less type; };

                template<class T>
                struct deref_t<T*>{ typedef double_deref_less type; };

                template<class T>
                void remove_unsorted_duplicates(std::list<T> &the_list) {
                  std::set<typename std::list<T>::iterator, typename deref_t<T>::type> found;
                  for (typename std::list<T>::iterator x = the_list.begin(); x != the_list.end();) {
                    if (!found.insert(x).second) {
                      x = the_list.erase(x);
                    }
                    else {
                      ++x;
                    }
                  }
                }

                struct EXTRACT_IF{
                    typedef std::list<infos_base*> result_type_single;
                    typedef std::list<infos_base*> result_type_all;
                    static void do_on_binary_trees(result_type_single & reslhs, result_type_single & resrhs,result_type_single & res){
                        res.merge(reslhs);
                        res.merge(resrhs);
                    }
                    static void do_on_unary_trees(result_type_single & ressub, result_type_single & res){
                        res.merge(ressub);
                    }
                    static void do_on_pred_true(infos_base* tree,result_type_single & res){
                        res.push_back(tree);
                    }
                    static void do_on_next_operation_merge(result_type_single & new_res, result_type_all & final_res){
                        final_res.merge(new_res);
                    }
                };


                template<class FILTER_T, class Pred>
                static typename FILTER_T::result_type_single filter(infos_base* const tree, Pred pred){
                    typedef typename FILTER_T::result_type_single res_t;
                    res_t res;
                    if(binary_tree_infos_base * p = dynamic_cast<binary_tree_infos_base *>(tree)){
                        res_t  reslhs(filter<FILTER_T,Pred>(&p->lhs(),pred));
                        res_t resrhs(filter<FILTER_T,Pred>(&p->rhs(),pred));
                        FILTER_T::do_on_binary_trees(reslhs,resrhs,res);
                    }
                    else if(unary_tree_infos_base * p = dynamic_cast<unary_tree_infos_base *>(tree)){
                        res_t ressub(filter<FILTER_T,Pred>(&p->sub(),pred));
                        FILTER_T::do_on_unary_trees(ressub,res);
                    }
                    if(pred(tree)){
                        FILTER_T::do_on_pred_true(tree,res);
                    }
                    return res;
                }

                template<class FILTER_T, class Pred>
                static typename FILTER_T::result_type_all filter(std::list<infos_base*> const & trees, Pred pred){
                    typedef typename FILTER_T::result_type_single res_t_single;
                    typedef typename FILTER_T::result_type_all res_t_all;
                    res_t_all res;
                    for(std::list<infos_base*>::const_iterator it = trees.begin() ; it != trees.end() ; ++it){
                        res_t_single tmp(filter<FILTER_T,Pred>(*it,pred));
                        FILTER_T::do_on_next_operation_merge(tmp,res);
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
                    return cast<T>(filter<EXTRACT_IF>(trees,is_type<T>));
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
