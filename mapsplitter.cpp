
#include <stdio.h>
#include <map>
#include <vector>

std::map<int, double> poses;
char tlf[1048576];
char blaj[1048576];
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Arguments should be mapfile coordfiles\n");
        return -1;
    }
    for (int i = 1; i < argc; i++)
    {
        FILE* in = fopen(argv[i], "rt");
        int chromo, pos;
        while (fscanf(in, "%d %d", &chromo, &pos) == 2)
        {
            poses[pos] = -1;
        }
    }
    int indcount = 0;

    int chromo, pos;
    double cm;

    while (scanf("%d %s %lf %d", &chromo, blaj, &cm, &pos) == 4)
    {
        auto iter = poses.find(pos);
        if (iter == poses.end())
        {
            fprintf(stderr, "Position %d unknown, skipping\n", pos);
        }

        iter->second = cm;
    }

    double cMperbp = 1e-6;
    int firstbp = -1;
    int lastbp = -1;
    double firstcM = -1;
    double lastcM = -1;
    double term;
    for (auto iter = poses.begin(); iter != poses.end(); iter++)
    {
        if (iter->second != -1)
        {
            firstbp = iter->first;
            firstcM = iter->second;
            break;
        }
    }

    for (auto iter = poses.begin(); iter != poses.end(); iter++)
    {
        if (iter->second != -1)
        {
            firstbp = iter->first;
            firstcM = iter->second;
        }
        if (lastbp <= firstbp)
        {
            auto jter = iter;
            for (jter++; jter != poses.end(); jter++)
            {
                if (jter->second != -1 && jter->first > firstbp)
                {
                    lastbp = jter->first;
                    lastcM = jter->second;
                    cMperbp = (lastcM - firstcM) / (lastbp - firstbp);
                    break;
                }
            }
        }

        if (iter->second == -1)
        {
            iter->second = firstcM + (iter->first - firstbp) * cMperbp;
        }
        
        if (iter == poses.begin())
        {
            term = iter->second;
        }
        iter->second -= term;
    }
    

    printf("%d\n", poses.size());
    for (auto pos : poses)
    {
        if (pos.second < 0)
        {
            fprintf(stderr, "No mapping position specified at %d\n", pos.first);
        }

        printf("%lf\n", pos.second);
    }
}