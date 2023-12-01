#include "main.h"

int main(int argc, char **argv)
{
    int numEle = 0;
    char *method = NULL, *problem = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-n", argv[index]))
        {
            numEle = atoi(argv[index + 1]);
        }
        if (strstr("-method", argv[index]))
        {
            method = argv[index + 1];
        }
        if (strstr("-problem", argv[index]))
        {
            problem = argv[index + 1];
        }
    }

    if (strstr(method, "trapezoid"))
    {
        if (strstr(problem, "task_1"))
        {
            TrapezoidIntegral(numEle, Task1Value);
        }
        else
        {
            TrapezoidIntegral(numEle, Task2Value);
        }
    }
    else if (strstr(method, "simpson"))
    {
        if (strstr(problem, "task_1"))
        {
            SimpsonIntegral(numEle, Task1Value);
        }
        else
        {
            SimpsonIntegral(numEle, Task2Value);
        }
    }

    return 0;
}

/*
 * command line
 * ./app_exe -n <num_1> -method <string_1> -problem <string_2>
 *     <num_1> number of elements
 *     <string_1> trapezoid or simpson
 *     <string_2> task_1 or task_2
 */