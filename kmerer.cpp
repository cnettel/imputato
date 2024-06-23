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
    int kmerlen;
    int count;

    if (argc < 3)
    {
        fprintf(stderr, "Expected arguments are KMERLEN COUNT");
        return -1;
    }

    if (sscanf(argv[1], "%d", &kmerlen) != 1)
    {
        fprintf(stderr, "Kmer len %s not an integer\n", argv[1]);
        return -1;
    }

    if (sscanf(argv[2], "%d", &count) != 1)
    {
        fprintf(stderr, "Count %s not an integer\n", argv[2]);
        return -1;
    }


    int pos;
    char tlf2[kmerlen + 1] = {0};

    while (scanf("%d %s", &pos, tlf) == 2)
    {        
        int len = strlen(tlf);
        int end = len - kmerlen;
        if (len - kmerlen <= 0)
        {
            fprintf(stderr, "Too short read %s at %d\n", tlf, pos);
            continue;
        }
        std::uniform_int_distribution start(0, end); 
        for (int k = 0; k < count; k++)
        {
            int thisstart = start(rng);
            strncpy(tlf2, &tlf[thisstart], kmerlen);
            printf("%d %s ", thisstart + pos, tlf2);
        }
        printf("\n");
    }
}
