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

template<class column> void dotransition(column& c, column& d, const map& themap, int marker, int d)
{

}

struct haplotype
{
    vector<bool> anyprior;
    vector<genprob> prior;
    vector<genprob> posterior;

    MatrixXf fwbw[2];
};

template<class column> void doemit(column& c, genprob& prior, int marker);
{

};

struct individ
{
    vector<int> genotypes;    
};