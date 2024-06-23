
#include <stdio.h>
#include <map>
#include <vector>

std::map<int, std::vector<int>> poses;
char tlf[1048576];
char blaj[1048576];
int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Arguments should be [1/2] coordfiles\n");
        return -1;
    }
    int copies = 0;
    if (sscanf(argv[1], "%d", &copies) != 1 || copies < 1 || copies > 2)
    {
        fprintf(stderr, "Copy count should be 1 or 2, not %s\n", copies);
        return -1;
    }
    for (int i = 2; i < argc; i++)
    {
        FILE* in = fopen(argv[i], "rt");
        int chromo, pos;
        while (fscanf(in, "%d %d", &chromo, &pos) == 2)
        {
            poses[pos] = {};
        }
    }
    int indcount = 0;

    while (fgets(tlf, 1048576, stdin))
    {
        if (tlf[0] == '#') continue;
        char attrs[1000];
        int chromo, pos;
        int strpos;
        int retval;
        if ((retval = sscanf(tlf, "%d %d %s %s %s %s %s %s GT %n", &chromo, &pos, blaj, blaj, blaj, blaj, blaj, blaj, &strpos)) != 8)
        {
            fprintf(stderr, "Read count %d\n", retval);
            continue;
        }

        int allele1, allele2, strpos2;
        auto iter = poses.find(pos);
        if (iter == poses.end())
        {
            fprintf(stderr, "Position %d unknown, skipping\n", pos);
        }
        else
        {
            fprintf(stderr, "Position %d known\n", pos);
        }
        while (sscanf(&tlf[strpos], "%d%c%d%n", &allele1, blaj, &allele2, &strpos2) == 3)
        {
            strpos += strpos2;
            if (blaj[0] != '|' && blaj[0] != ',' && blaj[0] != '/' && blaj[0] != '\\')
            {
                fprintf(stderr, "Unknown allele separator %c at %d\n", blaj[0], pos);
            }
            if (copies == 1)
            {
                iter->second.push_back(allele1 + allele2);
            }
            else
            {
                iter->second.push_back(allele1);
                iter->second.push_back(allele2);
            }                
        }
        if (indcount != 0 && iter->second.size() != indcount)
        {
            fprintf(stderr, "Inconsistent individual count %d vs. new %d\n", indcount, iter->second.size());
        }

        indcount = iter->second.size();
    }

    printf("%d\n", indcount);
    for (int i = 0; i < indcount; i++)
    {
        for (auto pos : poses)
        {
            int val = -1;
            if (i < pos.second.size())
            {
                val = pos.second[i];
            }
            printf("%d ", val);
        }
        printf("\n");
    }
}