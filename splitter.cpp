#include <stdio.h>
#include <string.h>
#include <vector>
#include <random>

using std::vector;
vector<char> data;
std::mt19937 rng;

char tlf[1048576];
int main(int argc, char** argv)
{
    int start;
    int end;
    int readlen;
    double coverage;
    if (argc < 5)
    {
        fprintf(stderr, "Expected arguments are STARTBP ENDBP READLEN COVERAGE");
        return -1;
    }

    if (sscanf(argv[1], "%d", &start) != 1)
    {
        fprintf(stderr, "Start position %s not an integer\n", argv[1]);
        return -1;
    }

    if (sscanf(argv[2], "%d", &end) != 1)
    {
        fprintf(stderr, "End position %s not an integer\n", argv[2]);
        return -1;
    }

    if (sscanf(argv[3], "%d", &readlen) != 1)
    {
        fprintf(stderr, "Read len %s not an integer\n", argv[3]);
        return -1;
    }

    if (sscanf(argv[4], "%lf", &coverage) != 1)
    {
        fprintf(stderr, "Read len %s not floating point value\n", argv[4]);
        return -1;
    }

    scanf("%s", tlf);    
    while (scanf("%s", tlf) == 1)
    {        
        for (char* tlf2 = tlf; *tlf2; tlf2++)
        {
            if (*tlf2 >= 'A' && *tlf2 <= 'z')
            {
                data.push_back(*tlf2);
            }
        }
    }

    double gap = readlen / coverage;
    std::poisson_distribution gaps(gap); 
    tlf[readlen] = 0;

    for (int now = start; now + readlen <= end; now += gaps(rng))
    {
        strncpy(tlf, &data[now], readlen);
        printf("%d %s\n", now, tlf);
        fprintf(stderr, "%d\n", now);
    }
    
}
