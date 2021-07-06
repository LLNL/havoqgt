// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <set>
#include <map>
#include <functional>
#include <iostream>
#include <sstream>
#include <nlohmann/json.hpp>

namespace clippy {
using integer = int64_t;
using number = double;
using boolean = bool;
using string = std::string;
using arraystring = std::vector<string>;
using arrayint = std::vector<int64_t>;
using arrayintint = std::vector<std::pair<int64_t, int64_t>>;
namespace detail {
template <typename T>
struct name_of_type : std::false_type {
  static constexpr const char *name = "UNKNOWN";
};

template <>
struct name_of_type<string> : std::true_type {
  static constexpr const char *name = "string";
};

template <>
struct name_of_type<integer> : std::true_type {
  static constexpr const char *name = "integer";
};

template <>
struct name_of_type<number> : std::true_type {
  static constexpr const char *name = "number";
};

template <>
struct name_of_type<boolean> : std::true_type {
  static constexpr const char *name = "bool";
};

template <>
struct name_of_type<arraystring> : std::true_type {
  static constexpr const char *name = "arraystring";
};

template <>
struct name_of_type<arrayint> : std::true_type {
  static constexpr const char *name = "arrayint";
};

template <>
struct name_of_type<arrayintint> : std::true_type {
  static constexpr const char *name = "arrayintint";
};

template <typename T>
inline std::string get_type_name() {
  // static_assert(
  //     name_of_type<T>::value,
  //     "Unsupported type, must be {int64_t, std::string, double, bool}");
  return std::string{name_of_type<T>::name};
}

}  // namespace detail

class clippy {
 public:
  clippy(const std::string &&name, const std::string &&desc) {
    m_json_config["method_name"] = name;
    m_json_config["desc"] = desc;
  }

  ~clippy() {
    if (return_values) { std::cout << m_json_return << std::endl; }
  }

  template <typename T>
  void add_required(const std::string &&name, const std::string &&desc) {
    add_required_validator<T>(name);
    size_t position = m_next_position++;
    m_json_config["args"][name]["type"] = detail::get_type_name<T>();
    m_json_config["args"][name]["desc"] = desc;
    m_json_config["args"][name]["position"] = position;
  }

  template <typename T>
  void add_optional(const std::string &&name, const std::string &&desc,
                    const T &default_val) {
    add_optional_validator<T>(name);
    m_json_config["args"][name]["type"] = detail::get_type_name<T>();
    m_json_config["args"][name]["desc"] = desc;
    m_json_config["args"][name]["position"] = -1;
    m_json_config["args"][name]["default_val"] = default_val;
  }

  template <typename T>
  void returns(const std::string &&desc) {
    m_json_config["returns"]["type"] = detail::get_type_name<T>();
    m_json_config["returns"]["desc"] = desc;
  }

  template <typename T>
  void to_return(const T &value) {
    // if (detail::get_type_name<T>() !=
    //     m_json_config["returns"]["type"].get<std::string>()) {
    //   throw std::runtime_error("clippy::to_return(value):  Invalid type.");
    // }
    return_values = true;
    m_json_return = value;
  }

  bool parse(int argc, char **argv) {
    const char *JSON_FLAG = "--clippy-help";
    const char *DRYRUN_FLAG = "--clippy-validate";
    if (argc == 2 && std::string(argv[1]) == JSON_FLAG) {
      std::cout << m_json_config;
      return true;
    }
    std::cin >> m_json_input;
    validate_json_input();

    if (argc == 2 && std::string(argv[1]) == DRYRUN_FLAG) { return true; }

    // Good to go for reals
    return false;
  }

  template <typename T>
  T get(const std::string &&name) {
    if (has_argument(name)) {  // if the argument exists
      if (detail::get_type_name<T>() !=
          m_json_config["args"][name]["type"].get<std::string>()) {
        throw std::runtime_error("clippy::get(name):  Invalid type.");
      }
      return m_json_input[name].get<T>();
    } else {  // it's an optional
              // std::cout << "optional argument found: " + name << std::endl;
      return m_json_config["args"][name]["default_val"].get<T>();
    }
  }

  bool has_argument(const std::string &name) {
    if (m_json_input.count(name) == 0) { return false; }
    return true;
  }

 private:
  void validate_json_input() {
    for (auto &kv : m_input_validators) { kv.second(m_json_input); }
    // TODO: Warn/Check for unknown args
  }

  template <typename T>
  void add_optional_validator(const std::string &name) {
    if (m_input_validators.count(name) > 0) {
      std::stringstream ss;
      ss << "CLIPPy ERROR:   Cannot have duplicate argument names: " << name
         << "\n";
      throw std::runtime_error(ss.str());
    }
    m_input_validators[name] = [name](const nlohmann::json &j) {
      if (j.count(name) == 0) { return; }  // Optional, only eval if present
      try {
        j[name].get<T>();
      } catch (const std::exception &e) {
        std::stringstream ss;
        ss << "CLIPPy ERROR:  Optional argument " << name << ": \"" << e.what()
           << "\"\n";
        throw std::runtime_error(ss.str());
      }
    };
  }

  template <typename T>
  void add_required_validator(const std::string &name) {
    if (m_input_validators.count(name) > 0) {
      throw std::runtime_error("Clippy:: Cannot have duplicate argument names");
    }
    m_input_validators[name] = [name](const nlohmann::json &j) {
      if (j.count(name) == 0) {
        std::stringstream ss;
        ss << "CLIPPy ERROR:  Required argument " << name << " missing.\n";
        throw std::runtime_error(ss.str());
      }
      try {
        j[name].get<T>();
      } catch (const std::exception &e) {
        std::stringstream ss;
        ss << "CLIPPy ERROR:  Required argument " << name << ": \"" << e.what()
           << "\"\n";
        throw std::runtime_error(ss.str());
      }
    };
  }

  nlohmann::json m_json_config;
  nlohmann::json m_json_input;
  nlohmann::json m_json_return;
  size_t m_next_position = 0;

  bool return_values = false;

  std::map<std::string, std::function<void(const nlohmann::json &)>>
      m_input_validators;
};
}  // namespace clippy
