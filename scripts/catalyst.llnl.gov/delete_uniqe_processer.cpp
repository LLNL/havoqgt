#include <iostream>
#include <cstdint>
#include <map>
#include <unordered_set>
#include <unordered_set>
#include <cassert>

using namespace std;

namespace std {
  template <> struct hash<std::pair<int64_t, int64_t>>
  {
    size_t operator()(const std::pair<int64_t, int64_t> & x) const
    {
        std::size_t h1 = std::hash<int64_t>()(x.first);
        std::size_t h2 = std::hash<int64_t>()(x.second);
        return h1 ^ (h2 << 1);
    }
  };
}

#if 1
int main() {
  unordered_set<std::pair<int64_t, int64_t>> table;

  while(1) {
    int64_t s;
    int64_t t;
    int64_t is_delete;

    cin >> s >> t >> is_delete;

    if (is_delete) {
      table.erase(std::pair<int64_t, int64_t>(s, t));
    } else {
      table.insert(std::pair<int64_t, int64_t>(s, t));
    }

    if ( cin.eof() ) { break; }
  }

  for ( const std::pair<int64_t, int64_t>& x: table ) {
    cout << x.first << " " << x.second << std::endl;
  }
  return 0;
}

#else

int main() {
  unordered_set<std::pair<int64_t, int64_t>> table;

  int64_t c_s = -1;

  while(1) {
    int64_t s;
    int64_t t;
    int64_t is_delete;

    cin >> s >> t >> is_delete;

    if (s != c_s) {
      s = c_s;

      for ( const std::pair<int64_t, int64_t>& x: table ) {
        cout << x.first << " " << x.second << std::endl;
      }
      table.clear();

      assert(!is_delete);
    }

    if (is_delete) {
      table.erase(std::pair<int64_t, int64_t>(s, t));
    } else {
      table.insert(std::pair<int64_t, int64_t>(s, t));
    }

    if ( cin.eof() ) { break; }
  }

  for ( const std::pair<int64_t, int64_t>& x: table ) {
    cout << x.first << " " << x.second << std::endl;
  }

  return 0;
}

#endif
