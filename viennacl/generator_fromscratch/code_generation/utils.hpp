#ifndef VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP
#define VIENNACL_GENERATOR_CODE_GENERATION_UTILS_HPP

#include <list>
#include <set>
#include <algorithm>
#include <ostream>
#include <map>

namespace viennacl{

    namespace generator{



        namespace code_generation{



            namespace utils{


                template<class T>
                struct is_type{
                    bool operator()(infos_base* p) const{
                        return dynamic_cast<T *>(p);
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

                class kernel_generation_stream : public std::ostream{
                private:
                    class kgenstream : public std::stringbuf{
                    public:
                        kgenstream(std::ostream& final_destination
                                   ,unsigned int const & tab_count) : final_destination_(final_destination)
                                                                      ,tab_count_(tab_count){ }
                        ~kgenstream() {  pubsync(); }
                        int sync() {
                            for(unsigned int i=0 ; i<tab_count_;++i)
                                final_destination_ << '\t';
                            final_destination_ << str();
                            str("");
                            return !final_destination_;
                        }
                    private:
                        std::ostream& final_destination_;
                        unsigned int const & tab_count_;
                    };

                public:
                    kernel_generation_stream(std::ostream& final_destination) : std::ostream(new kgenstream(final_destination,tab_count_))
                                                                                , tab_count_(0){ }
                    ~kernel_generation_stream(){ delete rdbuf(); }
                    std::string str(){
                        return static_cast<std::stringbuf*>(rdbuf())->str();
                    }

                    void inc_tab(){ ++tab_count_; }
                    void dec_tab(){ --tab_count_; }

                private:
                    unsigned int tab_count_;
                };

                struct EXTRACT_IF{
                    typedef std::list<infos_base*> result_type_single;
                    typedef std::list<infos_base*> result_type_all;
                    static void do_on_new_res(result_type_single & new_res, result_type_single & res){
                        res.merge(new_res);
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
                        FILTER_T::do_on_new_res(reslhs,res);
                        FILTER_T::do_on_new_res(resrhs,res);
                    }
                    else if(unary_tree_infos_base * p = dynamic_cast<unary_tree_infos_base *>(tree)){
                        res_t ressub(filter<FILTER_T,Pred>(&p->sub(),pred));
                        FILTER_T::do_on_new_res(ressub,res);
                    }
                    else if(function_base * p = dynamic_cast<function_base*>(tree)){
                        std::list<infos_base*> args = p->args();
                        for(std::list<infos_base*>::iterator it = args.begin() ; it!= args.end() ; ++it){
                            res_t newres(filter<FILTER_T,Pred>(*it,pred));
                            FILTER_T::do_on_new_res(newres,res);
                        }
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
                    return cast<T,infos_base>(filter<EXTRACT_IF>(trees,is_type<T>()));
                }


                template<class T>
                class cache_manager{
                public:
                    cache_manager( std::set<T *, viennacl::generator::deref_less>& expressions_read
                                  ,std::list<T *> const & expressions_write
                                  ,  utils::kernel_generation_stream & kss) : expressions_read_(expressions_read), expressions_write_(expressions_write)
                                                                              ,kss_(kss){
                    }

                    void fetch_entries(std::string const & idx){
                        for(typename std::set<T *, viennacl::generator::deref_less>::iterator it = expressions_read_.begin() ; it != expressions_read_.end() ; ++it){
                            T * p = *it;
                            p->access_name(p->name()+"_val");
                            kss_ << p->aligned_scalartype() << " " << p->generate() << " = " << p->name() << "[" << idx << "];" << std::endl;
                        }
                    }

                    void writeback_entries(std::string const & idx){
                        for(typename std::list<T * >::iterator it = expressions_write_.begin() ; it != expressions_write_.end() ; ++it){
                            T * p = *it;
                            kss_<< p->name() << "[" << idx << "]"<< " = "  << p->generate() << ";" << std::endl;
                        }
                    }

                private:
                    std::list<T * > expressions_write_;
                    utils::kernel_generation_stream & kss_;
                    std::set<T *, viennacl::generator::deref_less> & expressions_read_;

                };



            }

        }

    }

}
#endif // UTILS_HPP
