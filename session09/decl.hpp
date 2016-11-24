// handling symmetric cases

template <typename T, typename S>
struct Decl
{
    typedef typename Decl<S,T>::Type Type;
};

// case same types
template <typename T>
struct Decl<T, T>
{
    typedef T Type;
};

template <>
struct Decl<int, double>
{
    typedef double Type;
};
