// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_HDF5_VERTEX_DATA_WRITER_HPP
#define HAVOQGT_HDF5_VERTEX_DATA_WRITER_HPP

#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include <havoqgt/hdf5_utilities.hpp>

namespace havoqgt {

class hdf5_vertex_data_writer {
 public:
  template <typename vertex_id_type, typename vertex_data_type>
  using vertex_data_table_type = std::tuple<std::vector<vertex_id_type>,
                                            std::vector<vertex_data_type>>;

  explicit hdf5_vertex_data_writer(const std::string &hdf5_file_names)
      : m_ptr_file(nullptr) {
    priv_create_file(hdf5_file_names);
  }

  template <class vertex_id_type, class vertex_data_type>
  bool write_data(
      const vertex_data_table_type<vertex_id_type, vertex_data_type> &out_data,
      const std::string &vertex_column_name,
      const std::string &vertex_data_column_name) {
    auto success = priv_write_array_dataset(vertex_column_name,
                                            h5_data_type(vertex_id_type{}),
                                            std::get<0>(out_data));
    success &= priv_write_array_dataset(vertex_data_column_name,
                                        h5_data_type(vertex_data_type{}),
                                        std::get<1>(out_data));

    return success;
  }

 private:
  bool priv_create_file(const std::string &file_name) {
    const bool success = h5_try_and_catch([this, &file_name]() -> bool {
      m_ptr_file = std::make_unique<H5::H5File>(file_name.c_str(),
                                                H5F_ACC_TRUNC);
      return true;
    });
    if (!success) {
      m_ptr_file.reset(nullptr);
      return false;
    }
    return true;
  }

  template <typename T>
  bool priv_write_array_dataset(const std::string &   dataset_name,
                                const H5::DataType &  data_type,
                                const std::vector<T> &data) {
    if (!m_ptr_file) return false;

    return h5_try_and_catch([this, &dataset_name, &data_type,
                              &data]() -> bool {
      hsize_t       size = data.size();
      H5::DataSpace dataspace(1, &size);
      H5::DataSet   dataset = m_ptr_file->createDataSet(dataset_name.c_str(),
                                                        data_type,
                                                        dataspace);
      dataset.write(data.data(), data_type);
      return true;
    });
  }

  std::unique_ptr<H5::H5File> m_ptr_file;
};
}  // namespace havoqgt

#endif  // HAVOQGT_HDF5_VERTEX_DATA_WRITER_HPP
