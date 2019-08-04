#ifndef FUNCTION_TEST_H
#define FUNCTION_TEST_H
#include "functions.h"
#include "global.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <string>

TEST(UTILITY, SET_BIT)
{
    std::vector<uint64_t> flag;
    size_t idx = 0;
    ASSERT_DEATH(set_bit(idx, flag);, "");
    flag.push_back(0);
    idx = 1;
    set_bit(idx, flag);
    ASSERT_TRUE(get_bit(idx, flag));
    idx = 12;
    set_bit(idx, flag);
    ASSERT_TRUE(get_bit(idx, flag));
}
#endif // FUNCTION_TEST_H
