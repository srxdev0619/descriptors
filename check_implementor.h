#ifndef CHECK_IMPLEMENTOR_H
#define CHECK_IMPLEMENTOR_H
#include "descriptor34.h"
#include "checker.h"
#include <iostream>
#include "symmetric_functions.h"
#define max_check 2
#define epsilon_continuity 0.001


class check_implementor
{
private:
    checker *object_checker;
public:
    check_implementor();
    void check_all();
};

#endif // CHECK_IMPLEMENTOR_H
