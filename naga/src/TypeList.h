#ifndef OCTOPUS_TYPE_LIST_H
#define OCTOPUS_TYPE_LIST_H

#include <cstdint>
#include <iostream>

template <typename T>
static char const* name();

#define DECLARE_NAME(Type, Name)       \
  template <>                          \
  constexpr char const* name<Type>() { \
    return #Name;                      \
  }

DECLARE_NAME(int32_t, int32)
DECLARE_NAME(uint32_t, uint32)

DECLARE_NAME(float, float)
DECLARE_NAME(double, double)

// DECLARE_NAME(u32, Uint32)
// DECLARE_NAME(i32, Int32)

// DECLARE_NAME(u64, Uint64)
// DECLARE_NAME(i64, Int64)

// DECLARE_NAME(f16, Half);
// DECLARE_NAME(f32, Float);
// DECLARE_NAME(f64, Double);

template <typename T>
std::string appendName(std::string str) {
  return str + name<T>();
}

template <typename...>
struct TypeList;
using OCTypeList = TypeList<int, unsigned int, float,
                            double>;  //, u32, i32, u64, i64, f32, f64>;

template <typename Type, typename... Types>
struct TypeList<Type, Types...> {
  using Head = Type;
  using Tail = TypeList<Types...>;
};

template <>
struct TypeList<> {};

using EmptyList = TypeList<>;

template <typename Interfacer, typename TList>
struct TypeMap;

template <typename Interfacer>
struct TypeMap<Interfacer, EmptyList> {
  template <typename... Args>
  static void apply(Args&&...) {}
};

template <typename Interfacer, typename Type,
          typename... Types  //,
                             // typename TList = TypeList<Type, Types...>,
                             // typename Head = typename TList::Head, typename
                             // Tail = typename TList::Tail
          >
struct TypeMap<Interfacer, TypeList<Type, Types...>>
    : TypeMap<Interfacer, typename TypeList<Type, Types...>::Tail> {
  using TList = TypeList<Type, Types...>;
  using Head = typename TList::Head;
  using Tail = typename TList::Tail;

  using Base = TypeMap<Interfacer, Tail>;

  template <typename... Args>
  static void apply(Args&&... args) {
    Interfacer::template call<Head>(std::forward<Args>(args)...);
    Base::template apply(std::forward<Args>(args)...);
  }
};

template <typename Interfacer, typename TList, typename... Args>
void apply(Args&&... args) {
  TypeMap<Interfacer, TList>::template apply(std::forward<Args>(args)...);
}

// template< template <typename> class TemplateT, typename TypeList> struct
// TypeMap;

// template< template <typename> class TemplateT, typename... Types>
// struct TypeMap<TemplateT, TypeList<Types...>> : TemplateT<Types>...
// {};

#endif  // OCTOPUS_TYPE_LIST_H
