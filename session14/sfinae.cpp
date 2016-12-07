#include <cstdio>
#include <type_traits>


template <typename T>
struct Foo
{
};

template <typename Any>
struct IsFoo
{
    static constexpr bool value = false;
};

template <typename T>
struct IsFoo<Foo<T>>
{
    static constexpr bool value = true;
};

template <typename T>
struct Dummy
{
};

template <typename Any>
struct IsDummy {
    static constexpr bool value = false;
};

template <typename T>
struct IsDummy<Dummy<T>>
{
    static constexpr bool value = true;
};


template <typename T>
typename std::enable_if<IsFoo<T>::value,
         void>::type
print(const T &)
{
        std::printf("foo\n");
}

template <typename T>
typename std::enable_if<IsDummy<T>::value,
         void>::type
print(const T &)
{
        std::printf("dummy\n");
}

template <typename T>
typename std::enable_if<!IsDummy<T>::value && !IsFoo<T>::value,
         void>::type
print(const T &)
{
        std::printf("something else\n");
}


int
main()
{
    Foo<double>    foo;
    Dummy<float>   dummy;

    print(foo);
    print(dummy);
    print(10);
}
