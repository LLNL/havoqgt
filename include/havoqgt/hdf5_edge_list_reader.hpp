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

#include <functional>

#include "H5Cpp.h"

namespace havoqgt {

class hdf5_edge_list_reader {
 public:
  using vertex_id_type = int64_t;
  using weight_type    = int64_t;
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
    const bool success = priv_try_catch_h5([&file_nam, this]() -> bool {
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

  bool priv_read_array_dataset(const std::string &          dataset_name,
                               std::vector<vertex_id_type> *buf) {
    const bool success = priv_try_catch_h5([this, &dataset_name, buf]() {
      H5::DataSet dataset = m_ptr_file->openDataSet(dataset_name.c_str());

      hsize_t length = 0;
      dataset.getSpace().getSimpleExtentDims(&length);

      buf->resize(length);
      dataset.read(buf->data(), H5::PredType::NATIVE_INT64);

      return true;
    });

    return success;
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

  edge_list_type m_edge_buf;
  std::unique_ptr<H5::H5File> m_ptr_file;
};

}  // namespace havoqgt

#endif  // HAVOQGT_HDF5_EDGE_LIST_READER_INCLUDED
