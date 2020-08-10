#pragma once
namespace nucleona::range {
    template<class Rng>
    struct Iterator {
        using Type = std::decay_t<decltype(std::declval<Rng>().begin())>;
    };
    template<class Rng>
    using IteratorT = typename Iterator<Rng>::Type;


    template<class Rng>
    struct Reference {
        using Type = typename IteratorT<Rng>::reference;
    };
    template<class Rng>
    using ReferenceT = typename Reference<Rng>::Type;

    template<class Rng>
    struct Value {
        using Type = typename IteratorT<Rng>::value_type;
    };
    template<class Rng>
    using ValueT = typename Value<Rng>::Type;
}
