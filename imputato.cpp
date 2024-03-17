#include <algorithm>
#include <numeric>
#include <array>
#include <eigen3/Eigen/Dense>
#include <vector>

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

template<class column> void dotransition(column& c, column& d, const map& themap, int marker, int d);

struct haplotype
{
    vector<bool> anyprior;
    vector<genprob> prior;
    vector<genprob> posterior;

    MatrixXf fwbw[2];
};

vector<haplotype> haplotypes;

template<class column> void doemit(column& c, genprob& prior, int marker)
{
    for (int i = 0; i < haplotypes.size(); i++)
    {
        float val = 0.0f;
        for (int j = 0; j < 2; j++)
        {
            val += prior[j] * haplotypes[i][j];
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

const constexpr int ploidy = 4;

struct individ
{
    vector<int> genotypes;
    void samplehaplotypes(int index);
    void nudgehaplotypes(int index);
};

void individ::nudgehaplotypes(int index)
{
    for (int i = 0; i < genotypes.size(); i++)
    {
        array<float, ploidy> probs[2] = {{0.f}, {0.f}};
        probs[1][0] = 1.f;
        for (int j = 0; j < ploidy; j++)
        {
            std::fill(ploidy[j % 2].begin(), ploidy[j % 2].end(), 0.f);
            for (int k = 0; k <= j; k++)
            {
                float sum = haplotypes[index + j].posterior[i][0] + haplotypes[index + j].posterior[i][1];
                for (int l = 0; l < 2; l++)
                {
                    probs[!(j % 2)][k + l] += probs[!(j % 2)][k] * haplotypes[index + j].posterior[i][l];
                }
            }
        }
    }
}