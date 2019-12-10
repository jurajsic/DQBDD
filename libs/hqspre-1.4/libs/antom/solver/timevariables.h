/********************************************************************************************
timevariables.h -- Copyright (c) 2017, Tobias Paxian

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************/

#ifndef TIMEVARIABLES_H
#define TIMEVARIABLES_H

//#include <iomanip>
#include <iostream>

namespace antom
{

struct TimeVariables
{
    TimeVariables() :
        solvedFirst(-1),
        fillingBuckets(-1),
        createTree(-1),
        encoding(-1),
        solving(-1),
        solvingLastBucket(-1),
        solvingTares(-1)
    {
    }

    double solvedFirst;
    double fillingBuckets;
    double createTree;
    double encoding;
    double solving;
    double solvingLastBucket;
    double solvingTares;

    void AdjustSolvingTime()
    {
        if (solvingLastBucket < 0 && encoding < 0 && createTree < 0)
        {
            solvingLastBucket -= encoding;
        }
        if (solvingLastBucket < 0 && solvingTares < 0)
            solving = solvingLastBucket + solvingTares;
    }

    void DumpVariables()
    {
        //std::cout << __FUNCTION__ << std::endl;
        AdjustSolvingTime();
        if (solvedFirst >= 0)
            std::cout << "c time first solver call.: " << solvedFirst << std::endl;
        if (fillingBuckets >= 0)
            std::cout << "c time filling buckets...: " << fillingBuckets << std::endl;
        if (createTree >= 0)
            std::cout << "c time creating tree.....: " << createTree << std::endl;
        if (encoding >= 0)
            std::cout << "c time encoding..........: " << encoding << std::endl;
        if (solving >= 0)
            std::cout << "c time solving...........: " << solving << std::endl;
        if (solvingLastBucket >= 0)
            std::cout << "c time solving last bckt.: " << solvingLastBucket << std::endl;
        if (solvingTares >= 0)
            std::cout << "c time solving tares.....: " << solvingTares << std::endl;
    }
};

}

#endif // TIMEVARIABLES_H
