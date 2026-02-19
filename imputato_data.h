#ifndef IMPUTATO_DATA_H
#define IMPUTATO_DATA_H
#include <vector>
#include <array>
#include <eigen3/Eigen/Dense>

#ifndef IMPUTATO_SETTINGS
#include "imputato_settings.h"
#endif

using std::array;
using std::vector;

using genprob = array<float, 2>;
using Eigen::ArrayXXf;
using Eigen::ArrayXf;

using ratiotype = double;

struct map
{
    vector<unsigned int> chromstarts;
    vector<double> chromposes;
    vector<float> otherepses;
} ourmap;

struct haplotype
{
    vector<genprob> posterior;
    vector<genprob> posteriorwo;
    vector<float> sim;
    vector<array<float, ploidy>> crosssim;
    vector<float> desired;
    vector<float> offset;
    vector<float> momentum;
    vector<int> classes;
    array<array<float, numclasses>, numclasses> classweights;
    // should be weighted by number of haplotypes in "target" (outer), has to be symmetric
    // when excluding the weighting for a proper HMM
    // indexing is [to][from], for fast access in dotransition
    array<float, numclasses> initclassweights;

    array<int, 8>* allowedrefs = nullptr;

    ArrayXXf* fwbw;
    vector<double> renorm[2];
    double likelihood;
    genprob& getprior(int m) const;
    genprob& getnewprior(int m) const;
    float& getanyprior(int m) const;
    float& getnewanyprior(int m) const;
    int getindex() const;

    void dofwbw(bool fw, const map& themap, bool initatfw = true);
};

struct individ
{
    vector<int> genotypes;
    vector<array<int, 2>> reads;
    vector<array<genprob, ploidy>> rawprobs;
    array<int, 8> allowedrefs;
    vector<float> maxshared;
    vector<int> maxsharedid;
    array<float, ploidy + 1> genotypebias;
    bool initatfw;

    individ()
    {
        for (auto& allowed : allowedrefs)
        {
            allowed = -1;
        }
        for (auto& bias : genotypebias)
        {
            bias = 1;
        }
        initatfw = true;
    }
    void samplehaplotypes(int index);
    void nudgehaplotypes(int index);
    void doposteriorhaplotypes(int index);
    std::tuple<int, int, double> findflip(int index);
    bool handleflip(int index);
};

vector<individ> inds;

vector<vector<genprob> > priors;
vector<vector<genprob> > newpriors;
vector<vector<float> > anypriors;
vector<vector<float> > newanypriors;

vector<haplotype> haplotypes;

#endif