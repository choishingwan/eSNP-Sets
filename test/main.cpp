#include "global.h"
#include "gtest/gtest.h"
#include <iostream>
#include <string>
int main(int argc, char* argv[])
{
    // if (argc < 2) { std::cerr << "Usage: UnitTest <Data>" << std::endl; }
    // assert(argc >= 2);
    // path = argv[argc - 1];
    // just in case
    // path.append("/");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
