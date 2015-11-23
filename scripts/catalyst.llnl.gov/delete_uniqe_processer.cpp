#include <iostream>
#include <cstdint>
#include <map>
#include <unordered_set>

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

int main() {
#if 1
  unordered_set<std::pair<int64_t, int64_t>> table;
#else
  typedef multimap<int64_t,int64_t>::iterator itr_t;
  multimap<int64_t,int64_t> mlmap;
#endif

  while(1) {
    int64_t s;
    int64_t t;
    int64_t is_delete;
#if 1
    cin >> s >> t >> is_delete;
#else
    cin >> s >> t; 
    is_delete = false;
#endif
#if 1
    if (is_delete) {
      table.erase(std::pair<int64_t, int64_t>(s, t));
    } else {
      table.insert(std::pair<int64_t, int64_t>(s, t));
    }

#else
    if (is_delete) { // Delete operation
      itr_t itr = mlmap.find(s);
      for (itr_t end = mlmap.end(); itr != end; ++itr) {
       if (itr->first != s) break;
       if (itr->second == t) {
         mlmap.erase(itr);
       }
      }
    } else {  // Inserte operation
      bool is_unique = true;
      itr_t itr = mlmap.find(s);
      for (itr_t end = mlmap.end(); itr != end; ++itr) {
        if (itr->first != s) break;
        if (itr->second == t) {
          is_unique = false;
          break;
        }
      }
      if (is_unique)
        mlmap.insert(std::pair<int64_t, int64_t>(s, t));
    }
#endif

    //cout << s << "\t" << t << "\t" << is_delete << "\t" << endl;
    if ( cin.eof() ) { break; }
  }

#if 1
for ( const std::pair<int64_t, int64_t>& x: table ) {
  cout << x.first << "\t" << x.second << std::endl;
}
#else
//  itr_t itr = mlmap.begin();
//  for (itr_t end = mlmap.end(); itr != end; ++itr) {
//    cout << itr->first << "\t" << itr->second << std::endl;
//  }
#endif


  return 0;
}

