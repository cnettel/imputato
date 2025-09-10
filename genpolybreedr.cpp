#include <string>
#include <vector>
#include <map>

struct marker
{
    std::string name;
    int bp;
    double cM;
};

constexpr int ploidy = 4;
constexpr int readlimit = -1;
constexpr bool addrefsasunphased = false;
constexpr bool addrefsasphased = true;

bool fgetsvcf(char* line, int bytes, FILE* file)
{
    while (fgets(line, bytes, file))
    {
        if (line[0] != '#' || line[1] != '#') return true;
    }
    return false;
}

int main(int argc, char** argv)
{
    if (argc < 9)
    {
        fprintf(stderr, "Expected arguments pedfile.csv phasedfile.csv genfile.vcf chromonumber ref.hap inds.gen pop.map popgen.inds");
        return 1;
    }
    FILE* ped = fopen(argv[1], "r");
    FILE* phased = fopen(argv[2], "r");
    FILE* gen = fopen(argv[3], "r");
    if (!ped || !phased || !gen)
    {
        fprintf(stderr, "Unable to open some file\n");
        return 2;
    }
    int chromo;
    if (sscanf(argv[4], "%d", &chromo) != 1)
    {
        fprintf(stderr, "Invalid chromo number\n");
        return 3;
    }

    FILE* refout = fopen(argv[5], "w");
    FILE* genout = fopen(argv[6], "w");
    FILE* mapout = fopen(argv[7], "w");
    FILE* indsout = fopen(argv[8], "w");

    if (!refout || !genout || !mapout || !indsout)
    {
        fprintf(stderr, "Unable to open some file for writing\n");
        return 4;
    }

    std::vector<marker> markers;
    std::map<std::string, int> markermap; 

    // TODO INDS

    marker marker;
    char line[65536];
    char name[256];
    int pos;

    std::vector<std::string> pars;
    std::vector<std::vector<int> > refvals;


    fgets(line, 65536, phased);
    if (sscanf(line, "marker,chromosome,position,bp%n", &pos) < 4)
    {
        char* now = &line[pos];
        while (*now == ',')
        {
            now++;
            pos = 0;
            for (; *now != ',' && *now >= 32; name[pos++] = *now++)
            {}
            name[pos] = 0;
            pars.push_back(name);
            for (int i = 0; i < ploidy; i++)
            {
                refvals.emplace_back();
            }
        }
    }

    char chromostr[256];
    
    while (fgets(line, 16384, phased) && sscanf(line, "%[^,],%[^,],%lf,%d,%n", name, chromostr, &marker.cM, &marker.bp, &pos) >= 4)
    {
        int chromno;
        if (sscanf(chromostr, "chr%d", &chromno) != 1)
        {
            fprintf(stderr, "Unrecognized chromosome string %s\n", chromostr);
            continue;
        }

        if (chromno != chromo)
        {
            continue;
        }

        marker.name = name;
        markers.push_back(marker);
        markermap[marker.name] = markers.size() - 1;
        pos--;
        char* now = line + pos;

        for (int i = 0; i < pars.size(); i++, now += pos)
        {            
            int genos[ploidy];
            static_assert(ploidy == 4, "Wrong ploidy for hard-coded format string");
            if (sscanf(now, ",%d|%d|%d|%d%n", &genos[0], &genos[1], &genos[2], &genos[3], &pos) < ploidy)
            {
                fprintf(stderr, "Error reading parent genotype ind %s, marker %s:\n%s\n", pars[i].c_str(), marker.name.c_str(), now);
                return -1;
            }
            for (int j = 0; j < ploidy; j++)
            {
                refvals[i * 4 + j].push_back(genos[j]);
            }
        }
    }

    std::vector<std::string> inds;
    using genovector = std::vector<std::pair<int, std::pair<int, int>>>;
    std::map<std::string, genovector> indmap;
    fgets(line, 65536, ped);    
    while (fgets(line, 65536, ped))
    {
        int pop;
        int ploidy;
        char parnames[2][256];
        int res;
        if ((res = sscanf(line, "%[^,],%d,%[^,],%[^,],%d", name, &pop, parnames[0], parnames[1], &ploidy)) < 5)
        {
            fprintf(stderr, "Error reading pedigreee line (%d successful, maybe name %s):\n%s\n", res, res >= 1 ? name : "", line);
            continue;
        }

        if (pop < 1)
        {
            fprintf(stderr, "Skipping non-progeny individual %s\n", name);
        }
        inds.push_back(name);
        genovector& genos = indmap[name];
        genos.resize(markers.size(), {-1, {0, 0}});
    }

    fprintf(refout, "%d\n", refvals.size());
    for (auto& ref : refvals)
    {
        for (int i = 0; i < ref.size(); i++)
        {
            fprintf(refout, "%d%c", ref[i], i != ref.size() - 1 ? ' ' : '\n');
        }
    }

    int identified = 0;
    std::vector<genovector*> vcfgenos;
    fgetsvcf(line, 65536, gen);
    char* now = line;
    char str[256];
    for (int i = 0, pos; i < 9; i++, now += pos + 1)
    {       
        if (sscanf(now, "%[^\t\n\r]%n", str, &pos) < 1)
        {
            fprintf(stderr, "Error reading header column prefix %d\n", i);
            return -2;
        }
    }

    int identifiedinds = 0;
    for (int pos; sscanf(now, "%[^\t\n\r]%n", str, &pos) == 1 && pos; now += pos + 1)
    {
        auto iter = indmap.find(str);

        if (iter == indmap.end())
        {
            fprintf(stderr, "Unknown individual %s in geno vcf\n", str);
            vcfgenos.push_back(nullptr);
            continue;
        }
        vcfgenos.push_back(&iter->second);
        fprintf(indsout, "%s\n", str);
        identifiedinds++;
    }
    
    while (fgetsvcf(line, 65536, gen))    
    {
        int markerno = -1;
        now = line;
        for (int i = 0, pos; i < 9; i++,  now += pos + 1)
        {        
            if (sscanf(now, "%[^\t\n\r]%n", str, &pos) < 1)
            {
                fprintf(stderr, "Error reading vcf body column prefix %d\n", i);
                return -3;
            }

            if (i == 2)
            {
                auto iter = markermap.find(str);
                if (iter == markermap.end()) continue;
                markerno = iter->second;
                identified++;
            }            
        }

        if (markerno == -1) continue;

        for (int i = 0; i < vcfgenos.size(); i++, now += pos + 1)
        {
            if (sscanf(now, "%[^\t\n\r]%n", str, &pos) < 1)
            {
                fprintf(stderr, "Error reading vcf body ind column %d\n", i);
                return -4;
            }

            if (!vcfgenos[i]) continue;
            char tmp[255];
            int a, b;
            int dosage;
            int genos[ploidy];
            static_assert(ploidy == 4, "Wrong ploidy for hard-coded format string");
            if (sscanf(str, "%d/%d/%d/%d:%d,%d:", &genos[0], &genos[1], &genos[2], &genos[3], &a, &b) < 6)
            {
                fprintf(stderr, "Error reading genotype data at ind %d, marker %d:\n%s\n", i, markerno, str);
                continue;
            }
            dosage = genos[0] + genos[1] + genos[2] + genos[3];
            if (a + b < readlimit)
            {
                dosage = -1;
            }
            (*vcfgenos[i])[markerno] = {dosage, {a,b}};
        }
    }

    fprintf(genout, "%d\n", identifiedinds + (addrefsasunphased ? refvals.size() / 4 : 0) + (addrefsasphased ? refvals.size() : 0));
    for (genovector* gs : vcfgenos)
    {
        bool first = true;
        if (!gs) continue;

        for (auto g : *gs)
        {
            if (first)
            {
                first = false;
            }
            else
            {
                fprintf(genout, " ");
            }
            if (g.first == -1 && (g.second.first || g.second.second))
            {
                fprintf(genout, "%d;%d", g.second.first, g.second.second);
            }
            else
            {
                fprintf(genout, "%d", g.first);
            }
        }
        fprintf(genout, "\n");
    }

    if (addrefsasunphased)
    {
        for (int i = 0; i < refvals.size(); i += ploidy)
        {
            for (int m = 0; m < refvals[i].size(); m++)
            {
                int geno = 0;
                for (int j = 0; j < ploidy; j++)
                {
                    geno += refvals[i + j][m];
                }
                fprintf(genout, "%d%c", geno, m != refvals[i].size() - 1 ? ' ' : '\n');
            }
        }
    }

    if (addrefsasphased)
    {
        for (int i = 0; i < refvals.size(); i++)
        {      
            for (int m = 0; m < refvals[i].size(); m++)
            {      
                fprintf(genout, "%d%c", refvals[i][m] * ploidy, m != refvals[i].size() - 1 ? ' ' : '\n');
            }
        }
    }

    fprintf(mapout, "%d\n", markers.size());
    for (auto& m : markers)
    {
        fprintf(mapout, "%lf\n", m.cM);
    }

    printf("Parsed %d markers with %d refs, found %d in vcf for %d individuals\n", markers.size(), refvals.size(), identified, identifiedinds);
}