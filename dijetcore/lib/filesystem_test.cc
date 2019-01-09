#include "gtest/gtest.h"

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/string/string_utils.h"

#include <string>

#include "boost/filesystem.hpp"


using std::string;

TEST(CreateDir, TestSimple) {
  string directory_name = "TEST_EXAMPLE_DIRECTORY_3592852857";

  EXPECT_EQ(dijetcore::CreateDirectory(directory_name), true);

  boost::filesystem::path directory(directory_name);
  EXPECT_EQ(boost::filesystem::is_directory(directory), true);

  boost::filesystem::remove(directory);
  EXPECT_EQ(boost::filesystem::is_directory(directory), false);
}

TEST(CreateDir, TestTwoLevel) {
  string directory_name = "TEST_EXAMPLE_DIRECTORY_3592852857/INNER";
  string dir_inner = "INNER";
  string dir_outer = "TEST_EXAMPLE_DIRECTORY_3592852857";

  EXPECT_EQ(dijetcore::CreateDirectory(directory_name), true);

  boost::filesystem::path directory(directory_name);
  EXPECT_EQ(boost::filesystem::is_directory(directory), true);

  boost::filesystem::remove(directory);
  EXPECT_EQ(boost::filesystem::is_directory(directory), false);
  EXPECT_EQ(boost::filesystem::is_directory(dir_outer), true);

  boost::filesystem::path directory_2(dir_outer);
  boost::filesystem::remove(directory_2);
  EXPECT_EQ(boost::filesystem::is_directory(directory), false);
  EXPECT_EQ(boost::filesystem::is_directory(dir_outer), false);
}

TEST(Concatenation, Path) {
  string path1 = "first/path";
  string path2 = "second/path";
  string final_path = dijetcore::MakeString(path1, "/", path2);

  EXPECT_EQ(dijetcore::ConcatenatePath(path1, path2), final_path);
}

TEST(Concatenation, File) {
  string path1 = "first/path";
  string path2 = "file.cc";
  string final_path = dijetcore::MakeString(path1, "/", path2);

  EXPECT_EQ(dijetcore::ConcatenatePath(path1, path2), final_path);
}