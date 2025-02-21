#include <stdio.h>
#include <math.h>
#include <random>

std::mt19937 rng;

char tlf[1048576];
int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Specify a coverage and ploidy");
        return -1;
    }

    double coverage;
    int ploidy;
    if (sscanf(argv[1], "%lf", &coverage) != 1 || coverage < 0)
    {
        fprintf(stderr, "Argument %s not intepreted as a non-negative coverage.\n", argv[1]);
        return -1;
    }
    if (sscanf(argv[2], "%d", &ploidy) != 1 || ploidy < 1)
    {
        fprintf(stderr, "Argument %s not intepreted as a positive integral ploidy.\n", argv[2]);
        return -1;
    }

    int n;
    fgets(tlf, 1048576, stdin);
    if (sscanf(tlf, "%d", &n) != 1)
    {
        fprintf(stderr, "Unable to identify number of individuals from first line: %s", tlf);
        return -1;
    }

    printf("%d\n", n);

    for (int i = 0; i < n; i++)
    {
        fgets(tlf, 1048576, stdin);

        char* tlf2 = tlf;
        int tot = 0;
        int mismatches = 0;
        int pos;
        int geno;

        while (sscanf(tlf2, "%d%n", &geno, &pos) == 1)
        {
            tlf2 += pos;

            std::poisson_distribution p1(coverage / ploidy * geno);
            std::poisson_distribution p2(coverage / ploidy * (ploidy - geno));            
            printf("%d;%d ", p1(rng), p2(rng));
        }
        printf("\n");
    }
}