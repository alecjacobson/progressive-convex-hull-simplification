#include "test_framework.h"

#include <cstring>

// Runs every registered test. Optional argv[1] is a substring filter on the
// test name. Prints one line per test and a summary; exit code = #failures.
int main(int argc, char ** argv)
{
  const char * filter = (argc > 1) ? argv[1] : nullptr;

  int passed = 0, failed = 0, skipped = 0;
  for(const auto & tc : test_registry())
  {
    if(filter && tc.name.find(filter) == std::string::npos)
    {
      skipped++;
      continue;
    }
    try
    {
      tc.fn();
      printf("[ PASS ] %s\n", tc.name.c_str());
      passed++;
    }
    catch(const TestFailure & f)
    {
      printf("[ FAIL ] %s\n    %s\n", tc.name.c_str(), f.msg.c_str());
      failed++;
    }
    catch(const std::exception & e)
    {
      printf("[ FAIL ] %s\n    uncaught std::exception: %s\n",
             tc.name.c_str(), e.what());
      failed++;
    }
    catch(...)
    {
      printf("[ FAIL ] %s\n    uncaught unknown exception\n", tc.name.c_str());
      failed++;
    }
  }

  printf("\n%d passed, %d failed", passed, failed);
  if(skipped) printf(", %d skipped", skipped);
  printf("\n");
  return failed;
}
