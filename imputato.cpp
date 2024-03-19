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
        for (int m = 0; m < ploidy; m++)
        {
            array<float, ploidy> probs[2] = {{0.f}, {0.f}};
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
                        probs[!now][k + l] += probs[!now][k] * haplotypes[index + j].posterior[i][l] / sum;
                    }
                }

                
            }
        }

        float mean = 0.f;
        float sum = 1e-30f;
        for (int j = 0; j < ploidy; j++)
        {
            mean += j * probs[ploidy % 2][j];
            sum += probs[ploidy % 2][j];
        }

        mean /= sum;
    }
}