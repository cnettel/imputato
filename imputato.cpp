#include <algorithm>
#include <numeric>
#include <array>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <numeric>
#include <tuple>
#include <math.h>

using Eigen::ArrayXXf;
using Eigen::ArrayXf;
using std::array;
using std::vector;

using genprob = array<float, 2>;

struct map
{
    vector<int> chromstarts;
    vector<double> chromposes;  
} ourmap;

// Borrowed https://stackoverflow.com/questions/17719674/c11-fast-constexpr-integer-powers
constexpr int64_t ipow_(int base, int exp){
  return exp > 1 ? ipow_(base, (exp>>1) + (exp&1)) * ipow_(base, exp>>1) : base;
}
constexpr int64_t ipow(int base, int exp){
  return exp < 1 ? 1 : ipow_(base, exp);
}

const constexpr float Ne = 37.5;
const constexpr int ploidy = 4;
const constexpr int permcount = ipow(ploidy, ploidy);
float stepsize = 0.05;

template<class column> void doemit(column& c, genprob& prior, int marker);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d);

vector<vector<genprob> > priors;
vector<vector<genprob> > newpriors;
vector<vector<char> > anypriors;

struct haplotype
{
    vector<genprob> posterior;

    ArrayXXf* fwbw;
    vector<double> renorm[2];
    genprob& getprior(int m) const;
    genprob& getnewprior(int m) const;
    char& getanyprior(int m) const;
    int getindex() const;

    void dofwbw(bool fw, const map& themap)
    {
        ArrayXXf& myfwbw = fwbw[fw];
        int colcount = myfwbw.cols();

        int start = fw ? 0 : colcount - 1;
        int end = fw ? myfwbw.cols() : 0;
        int step = fw ? 1 : -1;
        int sidestep = fw ? 0 : -1;

        myfwbw.col(start).fill(1.0f / myfwbw.rows());
        renorm[fw][start] = 0.0f;

        for (int m = start; m != end; m += step)
        {
            auto col = myfwbw.col(m + sidestep);
            double srcrenorm = 0;
            if (m - sidestep - 1 >= 0)
            {
                int from = m - sidestep - 1;
                srcrenorm = renorm[fw][from];
                myfwbw.col(m + sidestep) = myfwbw.col(from);
                if (!fw && getanyprior(from)) doemit(col, getprior(from), from);
                dotransition(col, col, themap, from, step);
            }
            if (fw && getanyprior(m + sidestep)) doemit(col, getprior(m + sidestep), m + sidestep);

            for (int i = 0, j = getindex() / ploidy * ploidy; i < ploidy; i++, j++)
            {
                col(j) = 0;
            }

            float sum = col.sum();
            sum += 1e-30;
            
            renorm[fw][m + sidestep] = srcrenorm + log(sum);
            col *= expf(srcrenorm - renorm[fw][m + sidestep]);
        }
    }
};

vector<haplotype> haplotypes;

int haplotype::getindex() const
{
    return this - &haplotypes[0];
}

genprob& haplotype::getprior(int m) const
{
    return priors[m][getindex()];
}

genprob& haplotype::getnewprior(int m) const
{
    return newpriors[m][getindex()];
}

char& haplotype::getanyprior(int m) const
{
    return anypriors[m][getindex()];
}

int basehaps;

template<class column> void doemit(column& c, genprob& prior, int marker)
{
    vector<genprob>& ourPrior = priors[marker];
    vector<char>& ourAnyPrior = anypriors[marker];
    #pragma ivdep
    for (int i = 0; i < haplotypes.size(); i++)
    {
        float val = 0.0f;
        if (ourAnyPrior[i])
        {
            for (int j = 0; j < 2; j++)
            {
                val += prior[j] * ourPrior[i][j];
            }
        }
        else
        {
            /*for (int j = 0; j < 2; j++)
            {
                val += prior[j] * prior[j]; // Reasonable case is perfect match
            }*/
           val = 0.5;
        }
//        if (val < 0 || val > 1) printf("%f\n", val);
        c[i] *= val;
    }
}

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d)
{
    float dist = (themap.chromposes[marker + d] - themap.chromposes[marker]) * d * -0.02 * Ne;
    float nonrec = expf(dist);
    float rec = -expm1f(dist) / haplotypes.size();
    float sum = c.sum();

    c2 = c * nonrec + sum * rec;
}

std::mt19937 rng;

struct individ
{
    vector<int> genotypes;
    void samplehaplotypes(int index);
    void nudgehaplotypes(int index);
    void doposteriorhaplotypes(int index);
    std::tuple<int, int, double> findflip(int index);
    bool handleflip(int index);
};

vector<individ> inds;

void individ::samplehaplotypes(int index)
{
    // Very crude, biased
    std::uniform_real_distribution<float> distribution(0.99,1.01); 

    for (int j = 0; j < ploidy; j++)
    {        
        haplotypes[index + j].posterior.resize(genotypes.size());
        for (int i = 0; i < genotypes.size(); i++)
        {
            if (genotypes[i] >= 0)
            {
                float val = std::clamp((genotypes[i] / 1.0f / ploidy) * distribution(rng), 1e-5f, 1 - 1e-5f);
                haplotypes[index + j].getprior(i)[0] = 1.0f - val;
                haplotypes[index + j].getprior(i)[1] = val;
                haplotypes[index + j].getanyprior(i) = true;
            }
        }
    }
}

bool getploidyperm(int index, array<int, ploidy>& res)
{
    for (int k = 0; k < ploidy; k++)
    {
        res[k] = index % ploidy;
        index /= ploidy;

        for (int m = 0; m < k; m++)
        {
            if (res[m] == res[k]) return false;
        }
    }

    return true;
}

double likelihood;

std::tuple<int, int, double> individ::findflip(int index)
{
    int bestmarker = 0;
    int bestp = 0;
    double bestscore = -1.1e30f;

// TODO LESS MEMORY

    double scores[haplotypes[index].fwbw[0].cols()][permcount];
    double firstscore = 0;

    #pragma omp parallel for schedule(guided, 100), shared(firstscore)
    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        // TODO PRECALC ACCEL PLOIDY > 2
        bool first = true;
        double firstthisscore = 0;
        for (int p = permcount; p >= 0; p--)
        {
            array<int, ploidy> perm;
            if (!getploidyperm(p, perm))
            {
                scores[m][p] = -1.1e30f;
                continue;
            }

            /*if (first && m != 0)
            {
                scores[m][p] = -1.1e30f;
                first = false;
                continue;
            }*/            
            
            double sum = 0;
            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                //sum += haplotypes[index + j].renorm[1][m];
                sum += log((haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + perm[j]].fwbw[0].col(m)).sum() + 1e-30);
                //sum += haplotypes[index + perm[j]].renorm[0][m];
            }

            if (first)
            {
                if (m == 0)
                {
                    firstscore = sum;
                    for (int j = 0; j < ploidy; j++)
                    {
                        firstscore += haplotypes[index + j].renorm[1][m];   
                        firstscore += haplotypes[index + j].renorm[0][m];
                    }
                }
                sum += 0.01;
                firstthisscore = sum;
                first = false;
            }

            //if (index == 16) printf("Flip: %d %d %d %f\n", index, m, p, sum);
            scores[m][p] = sum - firstthisscore;
        }
    }

    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        for (int p = permcount; p >= 0; p--)
        {
            double sum = scores[m][p];
            if (sum < -1e30f) continue;

            if (sum > bestscore)
            {
                bestscore = sum;
                bestp = p;
                bestmarker = m;
            }
        }
    }    
    #pragma omp atomic
    likelihood += firstscore;

    return {bestmarker, bestp, bestscore};
}

bool individ::handleflip(int index)
{
    auto [bestmarker, bestp, bestscore] = findflip(index);


    array<int, ploidy> perm;
    bool straight = true;
    getploidyperm(bestp, perm);
    for (int j = 0; j < ploidy; j++)
    {
        straight &= perm[j] == j;
    }

    straight |= std::bernoulli_distribution()(rng);
    if (!straight) 
    {
        printf("Found flip for haplotype base index %d, marker %d, bestp %d, best score %f\n", index, bestmarker, bestp, bestscore);
        for (int j = 0; j < ploidy; j++)
        {
            printf("\t%d:%d", j, perm[j]);
        
        }
        printf("\n");
    }

    if (!straight)
    {
        array<genprob, ploidy> prior;
        array<genprob, ploidy> posterior;
        for (int i = 0; i < haplotypes[index].posterior.size(); i++)
        {
            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                prior[j] = haplotypes[index + j].getprior(i);
                posterior[j] = haplotypes[index + j].posterior[i];
            }

            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                int permval = i > bestmarker ? perm[j] : j;
                haplotypes[index + j].getnewprior(i) = prior[permval];
                haplotypes[index + j].posterior[i] = posterior[permval];
            }
        }
    }

    return !straight;
}

void individ::doposteriorhaplotypes(int index)
{
    ArrayXf probs;
#pragma omp parallel for private(probs), collapse(2), schedule(dynamic, 100)
    for (int j = 0; j < ploidy; j++)
    {
        for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
        {
            probs = haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + j].fwbw[0].col(m);
            haplotypes[index + j].posterior[m] = {0.0f, 0.0f};
            for (int k = 0; k < haplotypes.size(); k++)
            {
                //if (!haplotypes[k].getanyprior(k)) continue;

                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posterior[m][z] += priors[m][k][z] * probs(k) * anypriors[m][k];
                }
            }

            float sum = 1e-30f;
            for (int z = 0; z < 2; z++)
            {
                sum += haplotypes[index + j].posterior[m][z];
            }
            if (sum < 1e-10f) printf("HEJ %d %d %g\n", index + j, m, sum);
            sum = 1 / sum;
            for (int z = 0; z < 2; z++)
            {
                haplotypes[index + j].posterior[m][z] *= sum;
            }
        }
    }
}

void individ::nudgehaplotypes(int index)
{
#pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < genotypes.size(); i++)
    {
        if (genotypes[i] == -1) continue;

        int genotype = genotypes[i];

        {float probs[2] = {0.f};
        for (int m = 0; m < ploidy; m++)
        {
            for (int l = 0; l < 2; l++)
            {
                probs[l] += haplotypes[index + m].getprior(i)[l];
            }
        }

        float sum = probs[0] + probs[1] + 1e-30f;
        probs[0] /= sum;
        probs[1] /= sum;

        if (i < 10 && index < 4) printf(" %.3f/%.3f", probs[1]*ploidy, (float) genotype);}


        for (int m = 0; m < ploidy; m++)
        {
            array<float, ploidy + 1> probs[2] = {{0.f}, {0.f}};
            probs[1][0] = 1.f;
            int now = 1;
            for (int j = 0; j < ploidy; j++)
            {
                if (j == m)
                {
                    continue;
                }
                now = !now;
                std::fill(probs[now].begin(), probs[now].end(), 0.f);
                for (int k = 0; k < ploidy; k++)
                {
                    float sum = haplotypes[index + j].posterior[i][0] + haplotypes[index + j].posterior[i][1];
                    for (int l = 0; l < 2 && k + l < ploidy; l++)
                    {
                        probs[now][k + l] += probs[!now][k] * haplotypes[index + j].posterior[i][l] / sum;
                    }
                }
            }

            auto& priors = haplotypes[index + m].getprior(i);
            auto& newpriors = haplotypes[index + m].getnewprior(i);
            float val1 = (genotype ? probs[now][genotype - 1] : 0.f);
            float val2 = probs[now][genotype];
            //float diff = (val1 - val2) / (val1 + val2);

            float diff = logf((val1 + 1e-30f) / (val2 + 1e-30f));
            if (val1 + val2 < 1e-5f) printf("Stalemate at index %d, marker %d: %f, %f\n", index, i, val1, val2);
            for (int j = 0; j < 2; j++)
            {
                newpriors[j] = priors[j] * expf(diff * (j == 1 ? 1 : -1) * stepsize);
            }

            float sum = 0;
            for (int j = 0; j < 2; j++)
            {
                sum += newpriors[j];
            }

            sum = 1.f / sum;
            for (int j = 0; j < 2; j++)
            {
                newpriors[j] *= sum;
                newpriors[j] = std::clamp(newpriors[j], 1e-10f, 1.f);
            }
        }
    }
}

void initinds()
{
    int hapnum = haplotypes.size();
    basehaps = hapnum;
    haplotypes.resize(basehaps + inds.size() * ploidy);
    priors.resize(ourmap.chromposes.size());
    anypriors.resize(ourmap.chromposes.size());
    for (int i = 0; i < ourmap.chromposes.size(); i++)
    {
        priors[i].resize(haplotypes.size());
        anypriors[i].resize(haplotypes.size(), false);
    }

    for (individ& ind : inds)
    {
        ind.genotypes.resize(ourmap.chromposes.size());
        ind.samplehaplotypes(hapnum);
        hapnum += ploidy;
    }

    for (int h = basehaps; h < haplotypes.size(); h++)
    {
        auto& hap = haplotypes[h];
        for (int fw = 0; fw < 2; fw++)
        {
            hap.renorm[fw].resize(ourmap.chromposes.size());
        }        
    }

    newpriors = priors;
}

void doit()
{
    int hapnum = basehaps;
    ArrayXXf fwbw[ploidy][2];

    #pragma omp parallel for num_threads(3), private(hapnum, fwbw)
    for (int i = 0; i < inds.size(); i++)
    {
        #pragma omp parallel for num_threads(ploidy * 2), collapse(2), private(hapnum)
        for (int k = 0; k < ploidy; k++)
        {
            for (int fw = 0; fw < 2; fw++)
            {
                hapnum = basehaps + i * ploidy;
                haplotypes[hapnum + k].fwbw = fwbw[k];
                individ& ind = inds[i];
                haplotypes[hapnum + k].fwbw[fw].resize(haplotypes.size(), ourmap.chromposes.size());
                haplotypes[hapnum + k].dofwbw(fw, ourmap);
            }
        }
        hapnum = basehaps + i * ploidy;
        individ& ind = inds[i];
        bool flipped = ind.handleflip(hapnum);
        if (!flipped)
        {
            //printf("Nudge %d/%d\n", hapnum, haplotypes.size());
            ind.doposteriorhaplotypes(hapnum);
            ind.nudgehaplotypes(hapnum);
        }
    }

    priors = newpriors;
}

void readdummy(const char* mapname, const char* genoname)
{
    FILE* mapfile = fopen(mapname, "rt");
    ourmap.chromstarts.push_back(0);
    int d;
    fscanf(mapfile, "%d", &d);

    ourmap.chromposes.reserve(d);
    double prev = -1.f;
    for (int i = 0; i < d; i++)
    {
        double pos;
        fscanf(mapfile, "%lf", &pos);
        prev += 5e-5f;
        if (pos < prev)
        {
            pos = prev;
        }
        ourmap.chromposes.push_back(pos);
        prev = pos;
    }
    ourmap.chromstarts.push_back(d);

    FILE* indfile = fopen(genoname, "rt");
    int n;
    fscanf(indfile, "%d", &n);
    inds.resize(n);
    for (individ& ind : inds)
    {
        ind.genotypes.resize(d);
        for (int& g : ind.genotypes)
        {
            fscanf(indfile, "%d", &g);
        }
    }
}

void readrefs(const char* hapname)
{
    FILE* indfile = fopen(hapname, "rt");
    int n;
    fscanf(indfile, "%d", &n);
    priors.resize(ourmap.chromposes.size());
    anypriors.resize(ourmap.chromposes.size());
    for (int i = 0; i < ourmap.chromposes.size(); i++)
    {
        priors[i].resize(haplotypes.size() + n);
        anypriors[i].resize(haplotypes.size() + n, false);
    }
    
    for (int i = 0; i < n; i++)
    {
        haplotype& now = haplotypes.emplace_back();
        int d = ourmap.chromposes.size();
        now.posterior.resize(d);
        for (int j = 0; j < d; j++)
        {
            int val;
            fscanf(indfile, "%d", &val);
            if (val >= 0 && val <= 1)
            {
                now.getprior(j)[val] = 1.f - 1e-5f;
                now.getprior(j)[!val] = 1e-5f;
                now.getanyprior(j) = true;
            }
            else
            {
                now.getanyprior(j) = false;
            }
        }
    }
}

int main() 
{
    omp_set_max_active_levels(2);
    readdummy("potato_chr1.map", "potato_missing.gen");
    //inds.resize(2);
    initinds();
    for (int k = 0; k < 500; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 15; j++)
            {
                printf("%d %d", i, j);
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t%.3f %.3f", haplotypes[basehaps + i * ploidy + k].getprior(j)[1], haplotypes[basehaps + i * ploidy + k].getprior(j)[0]);
                }

                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t\t%.3f %.3f ", haplotypes[basehaps + i * ploidy + k].posterior[j][1], haplotypes[basehaps + i * ploidy + k].posterior[j][0]);
                }
                printf("\n");
            }
        }
        printf("Test! %d %lf\n", k, likelihood);
        likelihood = 0;
        doit();
    }

    FILE* out = fopen("potato.vcflike", "wt");
    for (int m = 0; m < ourmap.chromposes.size(); m++)
    {
        for (int i = 0; i < inds.size(); i++)
        {
            for (int k = 0; k < ploidy; k++)
            {
                //fprintf(out, "%c%.2f", k ? '|' : '\t', haplotypes[basehaps + i * ploidy + k].getprior(m)[1]);
                fprintf(out, "%c%.2f", k ? '|' : '\t', haplotypes[basehaps + i * ploidy + k].posterior[m][1]);
            }
        }
        fprintf(out, "\n");
    }
    fclose(out);

    out = fopen("potato.out", "wt");
    fprintf(out, "%d\n", inds.size());
    for (int i = 0; i < inds.size(); i++)
    {
        for (int m = 0; m < ourmap.chromposes.size(); m++)
        {
            int allele = 0;        
            for (int k = 0; k < ploidy; k++)
            {
                //fprintf(out, "%c%.2f", k ? '|' : '\t', haplotypes[basehaps + i * ploidy + k].getprior(m)[1]);
                allele += ((bool) (int) (haplotypes[basehaps + i * ploidy + k].posterior[m][1] * 2));
            }
            fprintf(out, "%d ", allele);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}
