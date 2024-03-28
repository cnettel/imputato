#include <algorithm>
#include <numeric>
#include <array>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>

using Eigen::MatrixXf;
using std::array;
using std::vector;

using genprob = array<float, 2>;

struct map
{
    vector<int> chromstarts;
    vector<float> chromposes;  
};

template<class column> void doemit(column& c, genprob& prior, int marker);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d);

struct haplotype
{
    vector<bool> anyprior;
    vector<genprob> prior;
    vector<genprob> posterior;

    // TODO: Fix major
    MatrixXf fwbw[2];

    void dofwbw(bool fw, const map& themap)
    {
        int colcount = fwbw[fw].cols();

        int start = fw ? 0 : colcount - 1;
        int end = fwbw[fw].cols();
        int step = fw ? 1 : -1;
        int sidestep = fw ? 0 : -1;

        fwbw[fw].col(start).fill(1.0f / fwbw[fw].rows());
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

void individ::doposteriorhaplotypes(int index)
{
    int bestmarker = 0;
    int bestp = 0;
    float bestscore = -std::numeric_limits<float>::infinity();
    for (int p = 0; p < permcount; p++)
    {
        array<int, ploidy> perm;
        if (!getploidyperm(p, perm)) continue;


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