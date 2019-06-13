#ifndef DIJETCORE_LIB_CONTAINERS_H
#define DIJETCORE_LIB_CONTAINERS_H

#include <limits>
#include <algorithm>

namespace dijetcore {

// finds the first occurence of object obj in container, and returns its index.
// returns std::numeric_limits<size_t>::max if obj is not in container
template <class Container, class T>
size_t FindFirst(const Container &container, const T &obj) {
  auto it = std::find(container.begin(), container.end(), obj);
  if (it == container.end())
    return std::numeric_limits<size_t>::max();
  else
    return std::distance(container.begin(), it);
}

} // namespace dijetcore

#endif // DIJETCORE_LIB_CONTAINERS_H
