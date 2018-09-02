#include "jet_analysis/util/arg_helper.hh"

#include <string>


using std::string;

// testing of arg_helper.hh functions, which are
// designed for command line argument parsing

int main() {
  
  // define a test string, test int, and test float
  string is_a_string = "hello";
  string is_an_int   = "12";
  string is_a_float  = "3.14";
  
  // test to see if CanCast() works as intended
  // (should be able to cast int as a float, but
  // do not allow floats to be shortened to ints)
  if (!CanCast<string>(is_a_string)) return 1;
  if (CanCast<double>(is_a_string))  return 1;
  if (!CanCast<int>(is_an_int))      return 1;
  if (!CanCast<double>(is_an_int))   return 1;
  if (!CanCast<double>(is_a_float))  return 1;
  if (CanCast<int>(is_a_float))      return 1;
  
  // change test strings to test Parse___Flag()
  is_a_string = "--isString=hello";
  is_an_int   = "--isInt=12";
  is_a_float  = "--isFloat=3.14";
  
  // define new strings to test ParseBoolFlag
  std::string is_false = "--testBool=false";
  std::string is_true = "--testBool=true";
  
  // initialize output
  string test_string = "";
  int test_int       = 0;
  double test_float  = 0.0;
  bool test_false    = true;
  bool test_true     = false;
  
  ParseBoolFlag(is_true, "--testBool", &test_true);
  ParseBoolFlag(is_false, "--testBool", &test_false);
  ParseStrFlag(is_a_string, "--isString", &test_string);
  ParseIntFlag(is_an_int, "--isInt", &test_int);
  ParseFloatFlag(is_a_float, "--isFloat", &test_float);
  
  // check for expected output from parse functions
  if (test_int != 12)         return 1;
  if (test_string != "hello") return 1;
  if (test_float != 3.14)     return 1;
  if (test_true != true)      return 1;
  if (test_false != false)    return 1;
  
  // test string parser
  string cast_to_double = "1.2,3424, 54.3";
  if (std::set<double>{1.2,3424,54.3} != ParseArgString<double>(cast_to_double)) return 1;
  
  string cast_to_int = "1, 534,45,   65";
  if (std::set<int>{1,534,45,65} != ParseArgString<int>(cast_to_int)) return 1;
  
  return 0;
}
