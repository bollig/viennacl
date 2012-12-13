#ifndef VIENNACL_TOOLS_ANY_HPP
#define VIENNACL_TOOLS_ANY_HPP

#include <typeinfo>

namespace viennacl{

    class any;

    template<class T>
    T any_cast(any& a);

    class value_base
    {
      public:
        virtual ~value_base() { }
        virtual value_base* clone() const = 0;
        virtual std::type_info const & type() const = 0;
    };

    template <class T>
    class value : public value_base
    {
        friend T any_cast<>(any& a);

        T t;

      public:
        value(const T& t_) : t(t_) { }
        value_base* clone() const
        {
            return new value(t);
        }

        std::type_info const &type() const {
            return typeid(T);
        }
    };

    class any
    {
        template<class T>
        friend T any_cast(any & a);

        value_base* v;

      public:
        any() : v(0) { }

        template <class value_type>
        any(const value_type& v_) : v(new value<value_type>(v_)) { }

        any(any const & other) : v(other.v ? other.v->clone() : 0) {      }

        any& operator=(const any& other)
        {
            if(&other != this)
            {
                any copy(other);
                swap(copy);
            }
            return *this;
        }

        void swap(any& other)
        {
            std::swap(v, other.v);
        }

        std::type_info const & type()
        {
          return v->type();
        }

        ~any() { delete v; }
    };

    class bad_any_cast : public std::bad_cast
    {
      public:
        virtual const char * what() const throw()
        {
            return "viennacl::bad_any_cast: "
                   "failed conversion using viennacl::any_cast";
        }
    };

    template <class T>
    T any_cast(any& a)
    {
      value<T>* v = dynamic_cast<value<T>*>(a.v);
      if(v == 0)
        throw bad_any_cast();
      else
        return v->t;
    }


}
#endif // ANY_HPP
