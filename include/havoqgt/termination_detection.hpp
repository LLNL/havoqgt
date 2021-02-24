// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_TERMINATION_DETECTION_HPP_INCLUDED
#define HAVOQGT_MPI_TERMINATION_DETECTION_HPP_INCLUDED

#include <havoqgt/mpi.hpp>
#include <utility>
#include <limits>

namespace havoqgt {

/*
Info about the states of termination detection -- this
is almost unreadable!

WAITING_FOR_PARENT_IRECV {  query_request & terminate_signal
if( local_received == local_retired) {
     if root, send msg to children
          else, MPI_Testany (query_request, terminate_signal)
}
-- locally al msgs have to be recived and sent
-- if root, ready
-- else, wait for parent signal
}

sent msg to children
WAITING_FOR_CHILDREN_ISEND  -- query_request & terminate_signal

WAITING_FOR_CHILDREN_IRECV{  -- query_response
  sum, (with local) and send to parent
    if root, eval and compare w/ 'previous' global reduce
}

WAIT_FOR_TERMINATE_CHILDREN_ISEND -- terminate_singal


irecv for parent (quiery and terminate)
irecv for childen (query response)
isend for children (query and terminate)
isend for parent (query response)

WAITING_FOR_INIT, WAITING_FOR_ISEND_CHILDREN, WAITING_FOR_RECV_CHILDREN, 
WAITING_FOR_ISEND_PARENT, 

*/

template<typename SizeType>
class termination_detection {
  private:
  typedef enum {WAITING_INIT, WAITING_ISEND_CHILDREN, WAITING_RECV_CHILDREN,
        WAITING_ISEND_PARENT} termination_detection_state_type;
  public:
  typedef SizeType size_type;
  typedef std::pair<size_type, size_type> status_response_type;

  termination_detection(MPI_Comm in_mpi_comm,
                                   int in_num_tree_children  = 2,
                                   int in_query_status_tag   = 2,
                                   int in_query_response_tag = 3,
                                   int in_terminate_tag      = 4) {
    m_mpi_comm = in_mpi_comm;
    m_num_tree_children = in_num_tree_children;
    m_query_status_tag = in_query_status_tag;
    m_query_response_tag = in_query_response_tag;
    m_terminate_tag = in_terminate_tag;
    CHK_MPI( MPI_Comm_rank( m_mpi_comm, &m_mpi_rank) );
    CHK_MPI( MPI_Comm_size( m_mpi_comm, &m_mpi_size) );
    m_num_waiting_recv_children = num_children();
    m_subtree_status_response = status_response_type(0,0);
    m_current_state = WAITING_INIT;
    m_previous_subtree_status_response = 
            status_response_type(std::numeric_limits<size_type>::max(),
                                 std::numeric_limits<size_type>::max());

    m_count_queued = 0;
    m_count_completed = 0;
    /*for(int i=0; i<m_mpi_size; ++i) {
      if(m_mpi_rank == i) {
    std::cout << m_mpi_rank << " parent_rank() = " << parent_rank() << std::endl;
    std::cout << m_mpi_rank << " begin_child_rank() = " << begin_child_rank() << std::endl;
    std::cout << m_mpi_rank << " is_leaf_rank() = " << is_leaf_rank() << std::endl;
    std::cout << m_mpi_rank << " num_children() = " << num_children() << std::endl;   
    std::cout << std::endl;
    }
    MPI_Barrier(in_mpi_comm);

    }
    exit(-1);*/
  }

  void inc_queued(size_t _i=1) { m_count_queued += _i; }
  void inc_completed(size_t _i=1) { m_count_completed += _i; }

  bool test_for_termination() {
    return test_for_termination_internal(m_count_queued, m_count_completed);
  }

  private:
  bool test_for_termination_internal(const size_type& in_queued,
                                     const size_type& in_completed) {
    switch(m_current_state) {
      case WAITING_INIT:
        return handle_waiting_init(in_queued, in_completed);
      case WAITING_ISEND_CHILDREN:
        return handle_waiting_isend_children(in_queued, in_completed);
      case WAITING_RECV_CHILDREN:
        return handle_waiting_recv_children(in_queued, in_completed);
      case WAITING_ISEND_PARENT:
        return handle_waiting_isend_parent(in_queued, in_completed);
    };
    return false;
  }

  bool handle_waiting_init(const size_type& in_queued,
                           const size_type& in_completed) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    m_subtree_status_response.first = 0;
    m_subtree_status_response.second = 0;
    if(m_mpi_rank == 0) {
      send_query_status_to_children();
      m_current_state = WAITING_ISEND_CHILDREN;      
    } else {
      if(mpi_iprobe(parent_rank(), m_terminate_tag)) {
        CHK_MPI( MPI_Recv(NULL, 0, MPI_BYTE, parent_rank(), m_terminate_tag, m_mpi_comm, MPI_STATUS_IGNORE) );
        send_terminate_to_children();
        m_current_state = WAITING_INIT;
        return true;
      }
      else if(mpi_iprobe(/*parent_rank()*/MPI_ANY_SOURCE, m_query_status_tag)) {
        //std::cout << m_mpi_rank << "ReceivedInit" << std::endl;
        CHK_MPI( MPI_Recv(NULL, 0, MPI_BYTE, parent_rank(), m_query_status_tag,
                          m_mpi_comm, MPI_STATUS_IGNORE) );
        send_query_status_to_children();
        m_current_state = WAITING_ISEND_CHILDREN;
      } 
    }
    return false;
  }

  bool handle_waiting_isend_children (const size_type& in_queued,
                                      const size_type& in_completed) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    while(!m_vec_req_isend_children.empty()) {
      if(mpi_test(m_vec_req_isend_children.back())) 
       m_vec_req_isend_children.pop_back();
      else 
        break;
    }
    if(m_vec_req_isend_children.empty()) 
      m_current_state = WAITING_RECV_CHILDREN;
    return false;
  }

  bool handle_waiting_recv_children (const size_type& in_queued,
                                     const size_type& in_completed) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    if(is_leaf_rank()) {
      m_subtree_status_response.first += in_queued;
      m_subtree_status_response.second += in_completed;
      isend_status_response_to_parent();
      m_current_state = WAITING_ISEND_PARENT;
    } else {
      status_response_type recv_buf;
      while(mpi_iprobe(MPI_ANY_SOURCE, m_query_response_tag)) {
        --m_num_waiting_recv_children;
        CHK_MPI( MPI_Recv( (void*) &recv_buf, sizeof(status_response_type), 
                           MPI_BYTE, MPI_ANY_SOURCE, m_query_response_tag, 
                           m_mpi_comm, MPI_STATUS_IGNORE) );
        m_subtree_status_response.first += recv_buf.first;
        m_subtree_status_response.second += recv_buf.second;
      }
      if(m_num_waiting_recv_children == 0) {
        m_num_waiting_recv_children = num_children();
        m_subtree_status_response.first += in_queued;
        m_subtree_status_response.second += in_completed;
        if(m_mpi_rank == 0) {
          if(m_subtree_status_response.first == m_subtree_status_response.second) {
            if(m_subtree_status_response == m_previous_subtree_status_response) {
              send_terminate_to_children();
              m_current_state = WAITING_INIT;
             return true;
            } else {
              m_previous_subtree_status_response = m_subtree_status_response;
            }
          } else {
            // This is a place we can debug termination detection
            //std::cout << "m_subtree_status_response = " << m_subtree_status_response.first
            //          << ", " << m_subtree_status_response.second << std::endl;
          }
          m_current_state = WAITING_INIT;
          return false;
        } else {
          isend_status_response_to_parent();
          m_current_state = WAITING_ISEND_PARENT;
        }
      }
    }
    return false;
  }

  bool handle_waiting_isend_parent (const size_type& in_queued,
                                    const size_type& in_completed) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    if(mpi_test(m_req_isend_parent)) {
      m_current_state = WAITING_INIT;
    }
    return false;
  }

  void send_terminate_to_children() {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    for(int i=0; i<num_children(); ++i) {
      int child_rank = i + begin_child_rank();
      CHK_MPI( MPI_Send(NULL, 0, MPI_BYTE, child_rank, m_terminate_tag,
                         m_mpi_comm) );
    }
  }

  bool mpi_test(MPI_Request& in_req) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    int flag(0);
    CHK_MPI( MPI_Test( &(in_req), &flag, MPI_STATUS_IGNORE) );
    return flag == 1;
  }

  void isend_status_response_to_parent() {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    CHK_MPI( MPI_Isend( (void*) &m_subtree_status_response, sizeof(status_response_type),
                        MPI_BYTE, parent_rank(), m_query_response_tag,
                        m_mpi_comm, &m_req_isend_parent) );
  }

  bool mpi_iprobe(int in_source, int in_tag) {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << " source = " << in_source << " tag = " << in_tag << std::endl;
    MPI_Status status;
    int flag(0);
    CHK_MPI( MPI_Iprobe(in_source, in_tag, m_mpi_comm, &flag, &status) );
    return flag == 1;
  }

  void send_query_status_to_children() {
    //std::cout << m_mpi_rank << " " << __FUNCTION__ << std::endl;
    for(int i=0; i<num_children(); ++i) {
      int child_rank = i + begin_child_rank();
      MPI_Request isend_request;
      //std::cout << m_mpi_rank << "MPI_Isend -- " << child_rank << " " << m_query_status_tag << std::endl;
      CHK_MPI( MPI_Isend(NULL, 0, MPI_BYTE, child_rank, m_query_status_tag,
                         m_mpi_comm, &isend_request) );
      m_vec_req_isend_children.push_back(isend_request);
    }
  }

  int parent_rank() { return (m_mpi_rank-1) / m_num_tree_children; }
  int begin_child_rank() { return (m_mpi_rank * m_num_tree_children)+1; }
  bool is_leaf_rank() { return begin_child_rank() >= m_mpi_size; }
  int num_children() {  
    int to_return = std::min(m_num_tree_children, m_mpi_size - begin_child_rank()); 
    if (to_return < 0) to_return = 0;
    return to_return;
  }


  ///Configuration parameters
  MPI_Comm m_mpi_comm;
  int m_num_tree_children;
  int m_query_status_tag;
  int m_query_response_tag;
  int m_terminate_tag;

  int m_mpi_size;
  int m_mpi_rank;

  termination_detection_state_type m_current_state;

  std::vector<MPI_Request> m_vec_req_isend_children;

  MPI_Request m_req_isend_parent;
  status_response_type* m_ptr_buf_isend_parent;
  int m_num_waiting_recv_children;
  status_response_type m_subtree_status_response;
  status_response_type m_previous_subtree_status_response;

  size_type m_count_queued;
  size_type m_count_completed;


  //MPI_Request m_
  //std::vector<MPI_Request> 

};

} //namespace havoqgt {

#endif //HAVOQGT_MPI_TERMINATION_DETECTION_HPP_INCLUDED
