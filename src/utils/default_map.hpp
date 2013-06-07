#pragma once
#include <map>
#include <iostream>

template<typename K, typename V> class DefaultMap {
public:
  typedef typename std::map<K, V>::iterator iterator;
  typedef typename std::map<K, V>::value_type value_type;
  typedef typename std::map<K,V>::const_iterator const_iterator;
  DefaultMap(const V& v): default_value(v) {}
  DefaultMap() {}
  V& operator[] (const K& k) {
    const_iterator it = m.find(k);
    if (it == m.end()) {
      m[k] = default_value;
    }
    return m[k];
  }
  V const& operator[] (const K& k) const {
    const_iterator it = m.find(k);
    if (it == m.end()) {
      return default_value;
    } else {
      return it->second;
    }
  }
  std::pair<iterator, bool> insert(const value_type& val) {
    return m.insert(val);
  }
  iterator begin() {
    return m.begin();
  }
  iterator end() {
    return m.end();
  }
  std::map<K, V> m;
  V default_value;
};

template class DefaultMap<std::string, double>;
