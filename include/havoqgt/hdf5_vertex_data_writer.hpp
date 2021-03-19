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

namespace havoqgt {

class hdf5_vertex_data_writer {
 public:
  using vertex_id_type   = uint64_t;
  using vertex_data_type = double;
  using vertex_data_table_type = std::tuple<std::vector<vertex_id_type>,
                                            std::vector<vertex_data_type>>;

  explicit hdf5_vertex_data_writer(const std::string &hdf5_file_names)
      : m_ptr_file(nullptr) {
    priv_create_file(hdf5_file_names);
  }

  bool write_data(const vertex_data_table_type &vertex_data) {
    auto success = priv_write_array_dataset(
        "vertex", H5::PredType::NATIVE_UINT64, std::get<0>(vertex_data));
    success &= priv_write_array_dataset("data", H5::PredType::IEEE_F64LE,
                                        std::get<1>(vertex_data));

    return success;
  }

 private:
  bool priv_create_file(const std::string &file_name) {
    const bool success = priv_try_catch_h5([this, &file_name]() -> bool {
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

    return priv_try_catch_h5([this, &dataset_name, &data_type,
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

  static bool priv_try_catch_h5(const std::function<bool()> &command,
                                const bool                   verbose = false) {
    if (!command) return true;

    try {
      H5::Exception::dontPrint();
      return command();
    } catch (const H5::FileIException &error) {
      if (verbose) error.printErrorStack();
      return false;
    } catch (const H5::DataSetIException &error) {
      if (verbose) error.printErrorStack();
      return false;
    } catch (const H5::DataSpaceIException &error) {
      if (verbose) error.printErrorStack();
      return false;
    } catch (const H5::DataTypeIException &error) {
      if (verbose) error.printErrorStack();
      return false;
    } catch (...) {
      if (verbose) std::cerr << "An exception was thrown" << std::endl;
      return false;
    }
    assert(false);
    return false;
  }

  std::unique_ptr<H5::H5File> m_ptr_file;
};
}  // namespace havoqgt

#endif  // HAVOQGT_HDF5_VERTEX_DATA_WRITER_HPP
