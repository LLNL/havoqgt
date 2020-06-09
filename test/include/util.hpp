#ifndef HAVOQGT_TEST_INCLUDE_UTIL_HPP
#define HAVOQGT_TEST_INCLUDE_UTIL_HPP

#include <string>
#include <cstdlib>

#define SUCCESS 0
#define FAILURE 1

namespace havoqgt::test {

/// \brief
/// The name of the environmental variable to set the data store directory
/// \Example
/// env HAVOQGT_TEST_DIR="/mnt/ssd/test" make test
const char *k_test_dir_env_name = "HAVOQGT_TEST_DIR";

/// \brief The default data store directory
const char *k_default_test_dir  = "/tmp";

namespace detail {
inline std::string get_test_dir_path() {
  if (const char *env_p = std::getenv(k_test_dir_env_name)) {
    return std::string(env_p);
  }
  return std::string(k_default_test_dir);
}
}  // namespace detail

inline std::string gen_test_dir_path(const std::string &dir_name) {
  return detail::get_test_dir_path() + "/" + dir_name;
}

}  // namespace havoqgt::test

#endif // HAVOQGT_TEST_INCLUDE_UTIL_HPP