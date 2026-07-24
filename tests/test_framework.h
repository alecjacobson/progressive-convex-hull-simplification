#pragma once
// Minimal zero-dependency test framework.
//
// Register a test with TEST(name) { ... } and assert with CHECK / CHECK_NEAR.
// run_tests.cpp provides main(): runs all registered tests (optionally filtered
// by a substring passed as argv[1]) and returns non-zero if any fail.
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

struct TestCase
{
  std::string name;
  std::function<void()> fn;
};

inline std::vector<TestCase> & test_registry()
{
  static std::vector<TestCase> r;
  return r;
}

struct TestRegistrar
{
  TestRegistrar(const std::string & name, std::function<void()> fn)
  {
    test_registry().push_back({name, std::move(fn)});
  }
};

struct TestFailure
{
  std::string msg;
};

#define TEST(name)                                                             \
  static void name();                                                          \
  static TestRegistrar test_reg_##name(#name, name);                           \
  static void name()

#define CHECK(cond)                                                            \
  do {                                                                         \
    if(!(cond))                                                                \
      throw TestFailure{std::string(__FILE__) + ":" +                          \
                        std::to_string(__LINE__) + "  CHECK failed: " #cond};  \
  } while(0)

#define CHECK_MSG(cond, msg)                                                   \
  do {                                                                         \
    if(!(cond))                                                                \
      throw TestFailure{std::string(__FILE__) + ":" +                          \
                        std::to_string(__LINE__) + "  " + (msg)};              \
  } while(0)

inline void check_near_impl(double a, double b, double tol, const char * ea,
                            const char * eb, const char * file, int line)
{
  if(std::isnan(a) || std::isnan(b) || std::fabs(a - b) > tol)
  {
    char buf[512];
    std::snprintf(buf, sizeof(buf),
      "%s:%d  CHECK_NEAR failed: %s = %.17g, %s = %.17g, |diff| = %.3g > %.3g",
      file, line, ea, a, eb, b, std::fabs(a - b), tol);
    throw TestFailure{buf};
  }
}

#define CHECK_NEAR(a, b, tol)                                                  \
  check_near_impl((a), (b), (tol), #a, #b, __FILE__, __LINE__)
