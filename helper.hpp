#ifndef HELPER_H
#define HELPER_H

#include <string>

using namespace std;

/* Helper function to check if a condition is met, similar to assert */
template <typename T>
T expect(T value, string error)
{
    if (!value)
        throw runtime_error(error);
    return value;
}

/* Helper function to cleanup a file */
void cleanup(const string &file)
{
    expect(!remove(file.c_str()), "Error deleting file " + file);
}

#endif