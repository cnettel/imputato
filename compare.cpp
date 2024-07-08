#include <stdio.h>
#include <math.h>

char tlf[2][1048576];
int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Specify two files to compare");
        return -1;
    }

    FILE* files[2];
    int n = 1 << 30;
    for (int i = 0; i < 2; i++)
    {
        if (!(files[i] = fopen(argv[i + 1], "rt")))
        {
            fprintf(stderr, "Unable to open file %d, %s\n", i + 1, argv[i + 1]);
            return -1;
        }

        fgets(tlf[i], 1048576, files[i]);
        int n2;
        sscanf(tlf[i], "%d", &n2);
        if (n2 < n) n = n2;
    }

    int totmismatch = 0;
    int totdisconc = 0;
    int totmissing = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            fgets(tlf[j], 1048576, files[j]);
        }

        char* tlf2[2] = {tlf[0], tlf[1]};
        int tot = 0;
        int mismatches = 0;
        int missing = 0;
        int disconc = 0;

        while (true)
        {
            int genos[2];
            bool ok = true;
            for (int i = 0; i < 2; i++)
            {
                int pos;
                if (sscanf(tlf2[i], "%d%n", &genos[i], &pos) != 1)
                {
                    if (tot < 10) fprintf(stderr, "Stopping at %d, %s\n", i, tlf[i]);
                    ok = false;
                }
                tlf2[i] += pos;
            }

            if (!ok) break;
            tot++;
            if (genos[0] == -1 || genos[1] == -1) missing++;
            else
            {
                mismatches += genos[0] != genos[1];
                disconc += abs(genos[0] - genos[1]);
            }
        }

        printf("Ind %d: %d/%d (%d missing, %d disconc)\n", i, mismatches, tot, missing, disconc);
        totmismatch += mismatches;
        totmissing += missing;
        totdisconc += disconc;
    }
    printf("Tot mismatches: %d\n", totmismatch);
    printf("Tot missing: %d\n", totmissing);
    printf("Tot disconc: %d\n", totdisconc);
}