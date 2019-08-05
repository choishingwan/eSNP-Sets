#ifndef FUNCTION_TEST_H
#define FUNCTION_TEST_H
#include "functions.h"
#include "global.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <random>
#include <string>

TEST(UTILITY, BIT_FLAG)
{
    std::vector<size_t> flag_empty;
    // check for empty vector
    // we always assume the vector to be initialized
    try
    {
        set_bit(0, flag_empty);
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    // check for size one vecgtor
    int idx;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<int> dist(-100, 100);
    std::uniform_real_distribution<size_t> size_dist(1, 1000);
    for (size_t j = 0; j < 50; ++j)
    {
        size_t cur_size = size_dist(mt);
        std::vector<size_t> flag(cur_size, 0);
        for (size_t i = 0; i < 50; ++i)
        {
            idx = dist(mt);
            if (idx < 0 || idx >= static_cast<int>((64 * cur_size)))

            { // when idx is less than 0, it will become a super large
                try
                {
                    set_bit(static_cast<size_t>(idx), flag);
                }
                catch (const std::out_of_range&)
                {
                    SUCCEED();
                }
            }
            else
            {
                set_bit(static_cast<size_t>(idx), flag);
                ASSERT_TRUE(get_bit(static_cast<size_t>(idx), flag));
            }
            // do a random set and get to check if it is valid
        }
    }
}

TEST(UTILITY, RANGE_BIT)
{
    std::vector<size_t> flag;
    // check for empty vector
    // we always assume the vector to be initialized
    try
    {
        range_set_bit(0, 1, flag);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    flag.push_back(0);
    try
    {
        range_set_bit(0, 65, flag);
        FAIL();
    }
    catch (const std::out_of_range&)
    {
        SUCCEED();
    }
    range_set_bit(0, 64, flag);
    for (size_t i = 0; i < 64; ++i) { ASSERT_TRUE(get_bit(i, flag)); }
}
TEST(UTILITY, NA_TEST)
{
    ASSERT_TRUE(is_na("NA"));
    ASSERT_TRUE(is_na("Na"));
    ASSERT_TRUE(is_na("NAN"));
    ASSERT_TRUE(is_na("NaN"));
    ASSERT_TRUE(is_na("na"));
    ASSERT_TRUE(is_na("nan"));
    ASSERT_TRUE(is_na("null"));
    ASSERT_TRUE(is_na("NULL"));
    ASSERT_FALSE(is_na("1NULL"));
    ASSERT_FALSE(is_na("123"));
    ASSERT_FALSE(is_na("-123"));
    ASSERT_FALSE(is_na("Null1"));
}
TEST(UTILITY, ATTRIBUTE)
{
    std::string gene_id, gene_name;
    bool res =
        parse_attribute("gene_id Test; gene_name Testing;", gene_id, gene_name);
    ASSERT_TRUE(res);
    ASSERT_STREQ(gene_id.c_str(), "Test");
    ASSERT_STREQ(gene_name.c_str(), "Testing");
    try
    {
        // malformed
        parse_attribute("gene_id; gene_name Testing;", gene_id, gene_name);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    gene_id = "";
    gene_name = "";
    res = parse_attribute("gene_id Test; ", gene_id, gene_name);
    ASSERT_FALSE(res);
    ASSERT_STREQ(gene_id.c_str(), "Test");
    ASSERT_STREQ(gene_name.c_str(), "");
    gene_id = "";
    gene_name = "";
    res = parse_attribute("gene_name Test; ", gene_id, gene_name);
    ASSERT_FALSE(res);
    ASSERT_STREQ(gene_name.c_str(), "Test");
    ASSERT_STREQ(gene_id.c_str(), "");
    gene_id = "";
    gene_name = "";
    res = parse_attribute("gene_id \"Test\"; gene_name \"Testing\";", gene_id,
                          gene_name);
    ASSERT_TRUE(res);
    ASSERT_STREQ(gene_id.c_str(), "Test");
    ASSERT_STREQ(gene_name.c_str(), "Testing");
}

TEST(UTILITY, GET_P)
{
    ASSERT_DOUBLE_EQ(get_p("0.1"), 0.1);
    ASSERT_DOUBLE_EQ(get_p("NA"), 2.0);
    ASSERT_DOUBLE_EQ(get_p("1000"), 1000);
    try
    {
        get_p("Test");
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
}

TEST(UTILITY, IS_GZ_FILE)
{
    GZSTREAM_NAMESPACE::ogzstream file;
    file.open("DEBUG.gz");
    file.close();
    ASSERT_TRUE(is_gz_file("DEBUG.gz"));
    std::remove("DEBUG.gz");
    file.open("DEBUG");
    file.close();
    ASSERT_TRUE(is_gz_file("DEBUG"));
    std::remove("DEBUG");
    try
    {
        is_gz_file("DEBUG");
        // not found
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    // should be robust against the suffix
    std::ofstream test;
    test.open("DEBUG.gz");
    test.close();
    ASSERT_FALSE(is_gz_file("DEBUG.gz"));
    std::remove("DEBUG.gz");
}
TEST(UTILITY, GET_BIM)
{
    ASSERT_STREQ(get_bim_name("Test").c_str(), "Test.bim");
    ASSERT_STREQ(get_bim_name("Test.bim").c_str(), "Test.bim");
    // in case someone want to name their plink file with bim
    ASSERT_STREQ(get_bim_name("bim").c_str(), "bim.bim");
    ASSERT_STREQ(get_bim_name("BIM").c_str(), "BIM.bim");
    // we assume the bim file is always well formated with a lower case
    // suffix.
    ASSERT_STREQ(get_bim_name("Test.BIM").c_str(), "Test.BIM.bim");
}

TEST(UTILITY, OPEN_GZ)
{
    GZSTREAM_NAMESPACE::igzstream gz;
    std::ifstream file;
    try
    {
        open_file("Test.gz", gz, file);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        // GZ file not found
        SUCCEED();
    }
    try
    {
        open_file("Test", gz, file);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        // normal file not found
        SUCCEED();
    }
    std::ofstream test;
    test.open("DEBUG");
    test.close();
    bool res = open_file("DEBUG", gz, file);
    // not gz file
    ASSERT_FALSE(res);
    ASSERT_TRUE(file.is_open());
    file.close();
    file.clear();
    gz.close();
    gz.clear();
    std::remove("DEBUG");
    GZSTREAM_NAMESPACE::ogzstream gz_test;
    gz_test.open("DEBUG.gz");
    gz_test.close();
    res = open_file("DEBUG.gz", gz, file);
    ASSERT_TRUE(res);
    ASSERT_TRUE(gz.good());
    ASSERT_FALSE(file.is_open());
    std::remove("DEBUG.gz");
    file.close();
    file.clear();
    gz.close();
    gz.clear();
    // test malform ending (should be robust against that)
    test.open("DEBUG.gz");
    test.close();
    res = open_file("DEBUG.gz", gz, file);
    ASSERT_FALSE(res);
    ASSERT_TRUE(file.is_open());
    std::remove("DEBUG.gz");
}
#endif // FUNCTION_TEST_H
