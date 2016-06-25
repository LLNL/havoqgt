#include <iostream>
#include <unistd.h>

int main(void)
{
  size_t sz;
  std::cin >> sz;
  size_t n = sz*1024*1024*1024;

  unsigned char* mem = new unsigned char[n];
  for (size_t i = 0; i < n; i+=4096) {
    mem[i] = 0;
  }

  std::cout << "Allocated(GB) " << sz << std::endl;


  while (1) {
    sleep(10000000);
  }

  return 0;
}