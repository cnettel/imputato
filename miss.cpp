#include <stdio.h>
#include <math.h>
#include <random>

std::mt19937 rng;

char tlf[1048576];
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Specify a missingness rate");
        return -1;
    }

    double missingness;
    if (sscanf(argv[1], "%lf", &missingness) != 1 || missingness < 0 || missingness > 1)
    {
        fprintf(stderr, "Argument %s not intepreted as a missingness rate from 0 to 1.\n", argv[1]);
        return -1;
    }
    std::bernoulli_distribution missing(missingness);

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

            if (missing(rng))
            {
                geno = -1;
            }
            printf("%d ", geno);
        }
        printf("\n");
    }
}