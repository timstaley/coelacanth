#include <UnitTest++/UnitTest++.h>
#include <iostream>
using namespace std;

//Global variable to tell the linked unit tests where to find resources:
string lucky_lib_test_resources_dir;


int main(int argc, char** argv)
{
    //comment
    cout << "Path to test resources: " << argv[1] <<endl;

    lucky_lib_test_resources_dir = string(argv[1]);

    lucky_lib_test_resources_dir+="/";

    return UnitTest::RunAllTests();
}

