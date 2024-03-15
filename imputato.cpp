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
        c2[i] = c[i] 
    }
}

struct individ
{
    vector<int> genotypes;    
};