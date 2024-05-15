#include <algorithm>
#include <numeric>
#include <array>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <numeric>
#include <tuple>

using Eigen::MatrixXf;
using std::array;
using std::vector;

using genprob = array<float, 2>;

struct map
{
    vector<int> chromstarts;
    vector<float> chromposes;  
} ourmap;

float Ne = 4000;

template<class column> void doemit(column& c, genprob& prior, int marker);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d);

struct haplotype
{
    vector<bool> anyprior;
    vector<genprob> prior;
    vector<genprob> posterior;

    // TODO: Fix major
    MatrixXf fwbw[2];
    vector<float> renorm[2];

    void dofwbw(bool fw, const map& themap)
    {
        int colcount = fwbw[fw].cols();

        int start = fw ? 0 : colcount - 1;
        int end = fwbw[fw].cols();
        int step = fw ? 1 : -1;
        int sidestep = fw ? 0 : -1;

        fwbw[fw].col(start).fill(1.0f / fwbw[fw].rows());
        renorm[fw][start] = 0.0f;
        for (int m = start; m != end; m += step)
        {
            auto col = fwbw[fw].col(m + sidestep);
            if (m)
            {
                int from = m - sidestep - 1;
                auto srccol = fwbw[fw].col(from);
                dotransition(srccol, col, themap, from, step);
            }
            doemit(col, prior[m + sidestep], m + sidestep);

            float sum = col.sum();
            sum += 1e-30;
            float oldrenorm = renorm[fw][start];
            renorm[fw][start] += log(sum);

            col *= exp(oldrenorm - renorm[fw][start]);
        }
    }
};

vector<haplotype> haplotypes;

template<class column> void doemit(column& c, genprob& prior, int marker)
{
    for (int i = 0; i < haplotypes.size(); i++)
    {
        float val = 0.0f;
        for (int j = 0; j < 2; j++)
        {
            val += prior[j] * haplotypes[i][marker][j];
        }
    }
}

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d)
{
    float dist = (map.chromposes[marker + d] - map.chromposes[marker]) * d;
    float nonrec = expf(dist * -0.02 * Ne);
    float rec = (1.0f - nonrec) / haplotypes.size();
    float sum = std::reduce(c.begin(), c.end());

    for (int i = 0; i < haplotypes.size(); i++)
    {
        c2[i] = c[i] * rec + (sum - c[i]) * nonrec;
    }
}

// Borrowed https://stackoverflow.com/questions/17719674/c11-fast-constexpr-integer-powers
constexpr int64_t ipow_(int base, int exp){
  return exp > 1 ? ipow_(base, (exp>>1) + (exp&1)) * ipow_(base, exp>>1) : base;
}
constexpr int64_t ipow(int base, int exp){
  return exp < 1 ? 1 : ipow_(base, exp);
}

const constexpr int ploidy = 4;
const constexpr int permcount = ipow(ploidy, ploidy);
std::mt19937 rng;

struct individ
{
    vector<int> genotypes;
    void samplehaplotypes(int index);
    void nudgehaplotypes(int index);
    void doposteriorhaplotypes(int index);
    std::tuple<int, int, float> findflip(int index);
    bool handleflip(int index);
};

void individ::samplehaplotypes(int index)
{
    // Very crude, biased
    std::uniform_real_distribution<float> distribution(0.99,1.01); 

    for (int j = 0; j < ploidy; j++)
    {
        haplotypes[index + j].posterior.resize(genotypes.size());
        haplotypes[index + j].prior.resize(genotypes.size());
        for (int i = 0; i < genotypes.size(); i++)
        {
            if (genotypes[i] >= 0)
        {
            float val = genotypes[i] / 1.0f * distribution(rng);
            haplotypes[index + j].prior[i][0] = 1.0f - val;
            haplotypes[index + j].prior[i][1] = val;
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

std::tuple<int, int, float> individ::findflip(int index)
{
    int bestmarker = 0;
    int bestp = 0;
    float bestscore = -std::numeric_limits<float>::infinity();
    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        for (int p = 0; p < permcount; p++)
        {
            array<int, ploidy> perm;
            if (!getploidyperm(p, perm)) continue;

            float sum = 0;
            for (int j = 0; j < ploidy; j++)
            {
                sum += haplotypes[index + j].renorm[1][m];
                sum += logf((haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + perm[j]].fwbw[0].col(m)).sum() + 1e-30f);
                sum += haplotypes[index + j].renorm[1][index + perm[j]];
            }

            if (sum > bestscore)
            {
                bestscore = sum;
                bestp = p;
                bestmarker = m;
            }
        }
    }    

    return {bestmarker, bestp, bestscore};
}

bool individ::handleflip(int index)
{
    auto [bestmarker, bestp, bestscore] = findflip(index);

    printf("Found flip for haplotype base index %d, marker %d, bestp %d, best score %f\n", index, bestmarker, bestp, bestscore);

    array<int, ploidy> perm;
    bool straight = true;
    for (int j = 0; j < ploidy; j++)
    {
        printf("\t%d:%d", j, perm[j]);
        straight &= perm[j] == j;
    }
    printf("\n");

    if (!straight)
    {
        array<genprob, ploidy> prior;
        array<genprob, ploidy> posterior;
        for (int i = bestmarker + 1; i < haplotypes[index].prior.size(); i++)
        {
            for (int j = 0; j < ploidy; j++)
            {
                prior[j] = haplotypes[index + j].prior[i];
                posterior[j] = haplotypes[index + j].posterior[i];
            }

            for (int j = 0; j < ploidy; j++)
            {
                haplotypes[index + j].prior[i] = prior[perm[j]];
                haplotypes[index + j].posterior[i] = posterior[perm[j]];
            }
        }
    }

    return !straight;
}

void individ::doposteriorhaplotypes(int index)
{
    MatrixXf probs;
    for (int j = 0; j < ploidy; j++)
    {
        for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
        {
            probs = haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + j].fwbw[0].col(m);
            haplotypes[index + j].posterior[m] = {0.0f, 0.0f};
            for (int k = 0; k < haplotypes[index].fwbw[0].rows(); k++)
            {
                if (!haplotypes[index + j].anyprior[k]) continue;

                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posterior[m][z] += haplotypes[k].prior[m][z] * probs(k);
                }
            }

            float sum = 1e-30f;
            for (int z = 0; z < 2; z++)
            {
                sum += haplotypes[index + j].posterior[m][z];
            }
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
    for (int i = 0; i < genotypes.size(); i++)
    {
        if (genotypes[i] == -1) continue;

        int genotype = genotypes[i];

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
                    for (int l = 0; l < 2; l++)
                    {
                        probs[now][k + l] += probs[!now][k] * haplotypes[index + j].posterior[i][l] / sum;
                    }
                }
            }

            float diff = (genotype ? probs[now][genotype - 1] : 0.f) - probs[now][genotype];
            auto& priors = haplotypes[index + m].prior[i];
            for (int j = 0; j < 2; j++)
            {
                priors[j] *= expf(diff * (j == 1 ? 1 : -1));
            }

            float sum = 0;
            for (int j = 0; j < 2; j++)
            {
                sum += priors[j];
            }

            sum = 1.f / sum;
            for (int j = 0; j < 2; j++)
            {
                priors[j] *= sum;
            }
        }
    }
}