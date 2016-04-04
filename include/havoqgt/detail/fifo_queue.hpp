#include <deque>

namespace havoqgt { namespace detail {

template <typename Visitor>
class fifo_queue {
protected:
  std::deque<Visitor> m_data;
public:
  fifo_queue() { }

  bool push( Visitor const& task) {
    m_data.push_back( task );
    return  true;
  }

  void pop() { m_data.pop_front(); }
  
  Visitor const& top() { return m_data.front(); }

  size_t size() const { return m_data.size(); }

  bool empty() const { return m_data.empty(); }

  void clear() { m_data.clear(); }
};
} /*end detail ns*/ }/*end havoqgt ns*/
