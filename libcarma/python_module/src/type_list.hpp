#ifndef CARMA_TYPE_LIST_H
#define CARMA_TYPE_LIST_H

#include <utility>
#include <carma.h>


template <typename T>
constexpr static char const* explicit_name();

#define DECLARE_NAME(Type, Name)       \
  template <>                          \
  constexpr char const* explicit_name<Type>() { \
    return #Name;                      \
  }

// DECLARE_NAME(u8, Uint8)
// DECLARE_NAME(i8, Int8)

// DECLARE_NAME(u16, Uint16)
// DECLARE_NAME(i16, Int16)

// DECLARE_NAME(u32, Uint32)
// DECLARE_NAME(i32, Int32)

// DECLARE_NAME(u64, Uint64)
// DECLARE_NAME(i64, Int64)

DECLARE_NAME(int, int);

#ifdef CAN_DO_HALF
DECLARE_NAME(half, half);
#endif
DECLARE_NAME(float, float);
DECLARE_NAME(double, double);

DECLARE_NAME(cuFloatComplex, float_complex);

template <typename T>
std::string appendName(std::string str) {
  return str + explicit_name<T>();
}

template <typename...>
struct GenericTypeList;

template <typename Type, typename... Types>
struct GenericTypeList<Type, Types...> {
  using Head = Type;
  using Tail = GenericTypeList<Types...>;
};

template <>
struct GenericTypeList<> {};

using EmptyList = GenericTypeList<>;

template <typename Interfacer, typename TList>
struct TypeMap;

template <typename Interfacer>
struct TypeMap<Interfacer, EmptyList> {
  template <typename... Args>
  static void apply(Args&&...) {}
};

template <typename Interfacer, typename Type,
          typename... Types  //,
                             // typename TList = GenericTypeList<Type,
                             // Types...>, typename Head = typename TList::Head,
                             // typename Tail = typename TList::Tail
          >
struct TypeMap<Interfacer, GenericTypeList<Type, Types...>>
    : TypeMap<Interfacer, typename GenericTypeList<Type, Types...>::Tail> {
  using TList = GenericTypeList<Type, Types...>;
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

// template< template <typename> class TemplateT, typename GenericTypeList>
// struct TypeMap;

// template< template <typename> class TemplateT, typename... Types>
// struct TypeMap<TemplateT, GenericTypeList<Types...>> : TemplateT<Types>...
// {};

#endif  // CARMA_TYPE_LIST_H
