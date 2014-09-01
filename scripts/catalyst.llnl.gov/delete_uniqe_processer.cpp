#include <iostream>
#include <map>
#include <cstdint>

using namespace std;


int main() {
  typedef multimap<int64_t,int64_t>::iterator itr_t;
  multimap<int64_t,int64_t> mlmap;

  while(1) {
    int64_t s;
    int64_t t;
    int64_t f;

    cin >> s >> t >> f;

    //cout << s << "\t" << t << "\t" << f << "\t" << endl;

    if (f) {
      itr_t itr = mlmap.find(s);
      for (itr_t end = mlmap.end(); itr != end; ++itr) {
       if (itr->first != s) break;
       if (itr->second == t) {
         mlmap.erase(itr);
       }
      }
    } else {
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
    if ( cin.eof() ) { break; }
  }

  itr_t itr = mlmap.begin();
  for (itr_t end = mlmap.end(); itr != end; ++itr) {
    cout << itr->first << "\t" << itr->second << std::endl;
  }

  return 0;
}
