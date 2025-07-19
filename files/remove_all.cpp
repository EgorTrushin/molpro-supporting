#include <experimental/filesystem>

namespace filesys = std::experimental::filesystem;

extern "C"
int remove_all(const char* file)
{
  return filesys::remove_all(file) ? 1 : 0;
}