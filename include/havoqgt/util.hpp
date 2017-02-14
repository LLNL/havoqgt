#include <chrono>
#include <sstream>
#include <string>
#include <vector>

/**
  * Returns elapsed (monotonic) time.
  */
double getElapsedTimeSecond(
  std::chrono::time_point<std::chrono::steady_clock> startTime,
  std::chrono::time_point<std::chrono::steady_clock> endTime) {
  return ((endTime - startTime).count()) *
    std::chrono::steady_clock::period::num /
    static_cast<double>(std::chrono::steady_clock::period::den);
}

template <typename TokenType>
std::vector<TokenType> split(const std::string& line, const char delim) {
  std::vector<TokenType> tokens;
  std::string token;
  std::istringstream iss(line);

  while(std::getline(iss, token, delim)) {    
    tokens.push_back(std::stoull(token));
  }
  
  return tokens;
} 

std::vector<std::string> split(const std::string& line, const char delim) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream iss(line);

  while(std::getline(iss, token, delim)) {
    tokens.push_back(token);
  }

  return tokens;
}

template <typename TokenType>
std::vector<TokenType> split_char(const std::string& line, const char delim) {
  std::vector<TokenType> tokens;
  std::string token;
  std::istringstream iss(line);

  while(std::getline(iss, token, delim)) {
    tokens.push_back(static_cast<TokenType>(token.at(0)));
  }

  return tokens;
}

