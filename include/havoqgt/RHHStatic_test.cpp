//
//  RHHStatic_test.cpp
//  static
//
//  Created by IwabuchiMlab on 2014/10/04.
//  Copyright (c) 2014年 IwabuchiMlab. All rights reserved.
//


/// How to use
/// $ g++ std=c++11 RHHStatic_test.cpp
/// $ ./a.out

#include <iostream>
#include <cstdint>
#include "RHHStatic.hpp"

void insertion_result(RHH::UpdateErrors err)
{
  switch (err) {
    case RHH::UpdateErrors::kSucceed:
      std::cout << "succeed\n";
      break;
    case RHH::UpdateErrors::kDuplicated:
      std::cout << "duplicated\n";
      break;
    case RHH::UpdateErrors::kReachingFUllCapacity:
      std::cout << "reaching full capacity\n";
      break;
    case RHH::UpdateErrors::kLongProbedistance:
      std::cout << "long probe distance\n";
      break;
  }
  
}

void erasion_result(bool err)
{
  if (err) {
    std::cout << "Succeed\n";
  } else {
    std::cout << "No elems\n";
  }
}

int main (void)
{
  RHH::RHHStatic<uint64_t, unsigned char, 256> *rhh = new RHH::RHHStatic<uint64_t, unsigned char, 256>();
  
  std::cout << "-------insert--------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    insertion_result(rhh->insert_uniquely(i, 0));
  }
  
  std::cout << "-------allocate chained RHHStatic--------\n";
  RHH::RHHStatic<uint64_t, unsigned char, 256> *rhh2 = new RHH::RHHStatic<uint64_t, unsigned char, 256>();
  rhh2->m_next_ = rhh;
  std::cout << "-------insert--------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    insertion_result(rhh2->insert_uniquely(i, 0));
  }
  std::cout << "--------erase-------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    erasion_result(rhh2->erase(i));
  }
  
  std::cout << "-------allocate chained RHHStatic--------\n";
  RHH::RHHStatic<uint64_t, unsigned char, 256> *rhh3 = new RHH::RHHStatic<uint64_t, unsigned char, 256>();
  rhh3->m_next_ = rhh2;
  std::cout << "------insert---------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    insertion_result(rhh3->insert_uniquely(i, 0));
  }
  std::cout << "--------erase-------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    erasion_result(rhh3->erase(i));
  }
  
  delete rhh;
  delete rhh2;
  delete rhh3;
  
  return 0;
}