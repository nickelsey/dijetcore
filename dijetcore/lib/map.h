#ifndef DIJETCORE_LIB_MAP_H
#define DIJETCORE_LIB_MAP_H

#include <unordered_map>

namespace dijetcore {
    // by default a dijetcore_map is an unordered map with no template
    // specialization by default
    template<class Key, class T, class Hash = std::hash<Key>>
    using dijetcore_map = std::unordered_map<Key, T, Hash>;
   
    // nick: I find using pairs and enumeration classes as keys can be
    // helpful sometimes, so I am implementing these before they are 
    // explicitly needed. 
    struct PairHash;
    template<class pairX, class pairY, class T, class Hash = PairHash>
    using pair_map = std::unordered_map<std::pair<pairX, pairY>, T, Hash>;

    struct EnumClassHash;
    template<class E, class T, class Hash = EnumClassHash>
    using enum_class_map = std::unordered_map<E, T, Hash>;

    // if using enum classes and pairs as keys, need to specify
    // the hash function
    struct EnumClassHash {
        template <typename T>
        std::size_t operator()(T t) const {
            return static_cast<std::size_t>(t);
        }    
    };

    struct PairHash {
        template <typename T, typename U>
        std::size_t operator()(const std::pair<T, U> &x) const {
            return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
        }
    };
} // namespace dijetcore

#endif // DIJETCORE_LIB_MAP_H
