// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_HDF5_EDGE_LIST_READER_INCLUDED
#define HAVOQGT_HDF5_EDGE_LIST_READER_INCLUDED

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

#include "H5Cpp.h"

#include <havoqgt/hdf5_utilities.hpp>

namespace havoqgt {

template <class _vertex_id_type = int64_t, class _weight_type = double>
class hdf5_edge_list_reader {
 public:
  using vertex_id_type = _vertex_id_type;
  using weight_type = _weight_type;

  using edge_list_type = std::tuple<std::vector<vertex_id_type>,
                                    std::vector<vertex_id_type>,
                                    std::vector<weight_type>>;

  hdf5_edge_list_reader(const std::string &hdf5_file_names)
      : m_ptr_file(nullptr) {
    priv_open_file(hdf5_file_names);
  }

  bool read(const std::string &col_name0, const std::string &col_name1,
            const std::string &weight_name) {
    if (!priv_read_array_dataset(col_name0, &std::get<0>(m_edge_buf)) ||
        !priv_read_array_dataset(col_name1, &std::get<1>(m_edge_buf))) {
      return false;
    }

    if (!weight_name.empty()) {
      if (!priv_read_array_dataset(weight_name, &std::get<2>(m_edge_buf))) {
        return false;
      }
    }

    return true;
  }

  const edge_list_type &edges() const { return m_edge_buf; }

 private:
  bool priv_open_file(const std::string &file_nam) {
    const bool success = h5_try_and_catch([&file_nam, this]() -> bool {
      m_ptr_file = std::make_unique<H5::H5File>(file_nam.c_str(),
                                          H5F_ACC_RDONLY);
      return true;
    });
    if (!success) {
      m_ptr_file.reset(nullptr);
      return false;
    }
    return true;
  }

  template <typename T>
  bool priv_read_array_dataset(const std::string &dataset_name,
                               std::vector<T> *   buf) {
    const bool success = h5_try_and_catch([this, &dataset_name, buf]() {
      H5::DataSet dataset = m_ptr_file->openDataSet(dataset_name.c_str());

      hsize_t length = 0;
      dataset.getSpace().getSimpleExtentDims(&length);

      buf->resize(length);
      dataset.read(buf->data(), h5_data_type(T{}));

      return true;
    });

    return success;
  }

  edge_list_type m_edge_buf;
  std::unique_ptr<H5::H5File> m_ptr_file;
};

}  // namespace havoqgt

#endif  // HAVOQGT_HDF5_EDGE_LIST_READER_INCLUDED
