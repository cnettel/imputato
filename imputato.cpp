#include <algorithm>
#include <numeric>
#include <array>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <numeric>
#include <tuple>
#include <math.h>

using Eigen::ArrayXXf;
using Eigen::ArrayXf;
using std::array;
using std::vector;

using ratiotype = double;

using genprob = array<float, 2>;

struct map
{
    vector<int> chromstarts;
    vector<double> chromposes;  
    vector<float> otherepses;
} ourmap;

// Borrowed https://stackoverflow.com/questions/17719674/c11-fast-constexpr-integer-powers
constexpr int64_t ipow_(int base, int exp){
  return exp > 1 ? ipow_(base, (exp>>1) + (exp&1)) * ipow_(base, exp>>1) : base;
}
constexpr int64_t ipow(int base, int exp){
  return exp < 1 ? 1 : ipow_(base, exp);
}

float Ne = 750;
const constexpr float Neend = 1.0;
const constexpr float Nedecay = 0.99;
const constexpr float Nestep = 0.0f;
const constexpr bool newNed = false;
const constexpr int ploidy = 4;
const constexpr int maxreads = 20; 
const constexpr int permcount = ipow(ploidy, ploidy);

using ArrayXPf = Eigen::Array<float, Eigen::Dynamic, ploidy>;

float stepsize = 0.015;
bool burnin = false;

template<class column> void doemit(column& c, float& anyprior, genprob& prior, int marker, int* indices);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d, int index);

vector<vector<genprob> > priors;
vector<vector<genprob> > newpriors;
vector<vector<float> > anypriors;
vector<vector<float> > newanypriors;

constexpr bool disableplacement = true;
constexpr bool disableperm = true;
constexpr bool permpostburnin = false;
constexpr bool altperm = false;
constexpr bool mulpriorrest = true;
constexpr bool mulpriorself = false;
constexpr bool antipriorself = false;
constexpr bool antipriorpp = false;
//constexpr bool antidiffp = false;
constexpr bool propriorself = false;
constexpr bool allhets = false;
constexpr bool burninassgn = false;
constexpr bool postassgn = false;
constexpr bool earlypost = true;
constexpr bool allnotme = false;
constexpr bool selfposteriorwo = true;
constexpr bool restposteriorwo = true;
constexpr float momspeed = 0.02f;
constexpr float killflipm = 0.05f;
constexpr bool selfpriorunass = false;
constexpr bool nonwoearly = false;
constexpr bool fullwo = true;
constexpr bool liftmeannprior = false;
constexpr bool antiselfmean = false;
double posteriormix = 0.0;
constexpr double respostmix = 0.0001;
double tension = 1.00;
constexpr double tensiongrow = 1.0000;
constexpr int tensionreset = 2000;
constexpr bool extremetension = false;
constexpr bool nocentertaper = false;
constexpr double offsetdecay = 1.000;
constexpr bool fillinmissingwo = true;
constexpr bool neutralmissing = false;
constexpr bool noanypriorweight = true;
constexpr bool updallpriors = false;
constexpr bool simplestep = false;
constexpr bool modifiedclamp = false;
constexpr bool domaxcentered = true;
constexpr bool tensionoffset = false;
constexpr bool nozeroone = false;
constexpr bool simoffset = true;
constexpr bool crosssimoffset = true;
constexpr double csscale = 10;
constexpr bool priororig = false;
constexpr bool priorpowo = true; // powoorig makes more sense
constexpr bool mixorig  = false;
constexpr bool postorig = false;
constexpr bool expklsim = false;
constexpr bool extremsim = false;
constexpr bool invcs = false;
constexpr bool simpowo = false;
constexpr double midpointcap = 30;
constexpr bool caponpriors = true;
constexpr bool simposteriormix = true;
constexpr bool antisimposterior = true;
constexpr bool binomsimcomp = false;
constexpr bool neverflip = false;
constexpr bool preextremis = false;
constexpr bool highsimpost = false;
constexpr bool antihisim = true;
constexpr bool postmixred = false;
constexpr bool postmixreset = false;
constexpr bool newpostmix = true;
constexpr double csbump = 0.000;
constexpr int endstepgrow = 1650;
constexpr int startstepshrink = 20000;
constexpr bool redcertmix = false;
constexpr bool simredcert = false;
constexpr bool antiredcert = true;
constexpr double certfactor = 4;
constexpr double certterm = 0;
constexpr bool killibd2trans = false;
constexpr bool killibd2unc = false;
constexpr bool genounc = false;
constexpr bool sepgenounc = false;
constexpr bool snowball = false;
constexpr bool bigsnow = true;
constexpr float momspeed2 = 0.99f;
constexpr bool clampmom = true;
constexpr double clamplim = 2.0;
constexpr bool arimeanmix = false;
constexpr bool logitmeanmix = true;
constexpr bool clampmix = false;
constexpr bool agnosticflip = false;
constexpr bool bothagn = false;
constexpr bool advstep = true;
constexpr double minstep = 1e-32;
constexpr bool selflipcs = false;
constexpr bool midpoint11m = false;
constexpr bool rngflwght = true;
double flipdrag = 0;
constexpr double fldrstep = 0.001;
constexpr bool plaincsdiff = false;
constexpr bool momstepclamp = false;
constexpr bool noburninmom = false;
constexpr int advpostsum = 1;
constexpr bool minimalstep = true;
constexpr bool guardminimum = true;
constexpr bool logitstep = true;
constexpr float dampextreme = 0.999f;
constexpr bool stepszoffs = true;
constexpr bool ibd2sort = true;
constexpr bool nodropysum = true;
constexpr int ibdgroupsize = 4;
constexpr bool sortemit = false;
constexpr bool resetsort = true;
constexpr bool sortsim = true;
constexpr bool sortflip = true;
constexpr bool progribd = true;
constexpr double levelibd = 3;
constexpr bool ibddecl = true;
constexpr bool adjibd = true;
constexpr double declstep = 1;
constexpr double minlevelibd = 3;
constexpr bool preindex = true;
constexpr double ibdfloor = 0;
constexpr bool orgsumw = true;
constexpr double ibdeps = 1e-9;
constexpr bool readj = true;
constexpr bool tightreadj = true;
constexpr bool dofixed = false;
constexpr bool groupibd2 = true;
constexpr bool antigroup = true;
constexpr bool antigflip = true;
constexpr float flipscale = 0.99f;
constexpr bool randflippos = false;
constexpr bool noiseflippos = false;
constexpr bool ibdfactors = true;
constexpr bool ibdmax = false;
constexpr bool ibdmax2 = false;
constexpr bool groupmax = false;
constexpr bool antifactors = true;
constexpr bool powofactors = false;
constexpr bool vetoflip = true;
constexpr bool crossibd = true;
constexpr bool nonsimfactor = true;
constexpr bool mulminibd = true;
constexpr bool dovar = true;
constexpr bool dovarunc = false;
constexpr bool scaleunc = false;
constexpr bool onlyvar = false;
constexpr double uncshift = 0.0;
constexpr bool scaleanyprior = false;
constexpr bool wounc = false;
constexpr bool scaleanypriorw = true;
constexpr bool nocshz = false;
constexpr bool onlyref = true;
constexpr float epsothergeno = 0 * 0.005 / (ploidy - 1);
constexpr bool domarkeps = true;
constexpr bool uncnogeno = true;
constexpr float epsiloncM = 5e-3f;
constexpr bool halfpar = false;
constexpr bool fixatone = false;
constexpr float offsmagn = 0.1;
constexpr float filterlevel = 1.0f;
constexpr float refeps = 1e-10f;
constexpr bool weakeneps = true;
constexpr bool rewcsdiff = false;
constexpr float updeps = 1e-5f;
constexpr bool filterrefs = true;
constexpr bool mulsimoffset = true;
constexpr bool halfparinit = true;

struct haplotype
{
    vector<genprob> posterior;
    vector<genprob> posteriorwo;
    vector<float> sim;
    vector<array<float, ploidy>> crosssim;
    vector<float> desired;
    vector<float> offset;
    vector<float> momentum;
    array<int, 8>* allowedrefs = nullptr;

    ArrayXXf* fwbw;
    vector<double> renorm[2];
    genprob& getprior(int m) const;
    genprob& getnewprior(int m) const;
    float& getanyprior(int m) const;
    float& getnewanyprior(int m) const;
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

        int indices[myfwbw.rows() / 2 / ibdgroupsize]; //ibd2sort

        for (int m = start; m != end; m += step)
        {
            auto col = myfwbw.col(m + sidestep);
            double srcrenorm = 0;
            if (m - sidestep - 1 >= 0)
            {
                int from = m - sidestep - 1;
                srcrenorm = renorm[fw][from];
                myfwbw.col(m + sidestep) = myfwbw.col(from);
                if (!fw /*&& getanyprior(from)*/) doemit(col, getanyprior(from), getprior(from), from, indices);
                dotransition(col, col, themap, from, step, getindex());
            }

            if (fw /*&& getanyprior(m + sidestep)*/)
            {
                if (fullwo)
                {
                    fwbw[2].col(m + sidestep) = myfwbw.col(m + sidestep);
                }
                doemit(col, getanyprior(m + sidestep), getprior(m + sidestep), m + sidestep, indices);
            }

            for (int i = 0, j = getindex() / ploidy * ploidy; i < ploidy; i++, j++)
            {
                col(j * 2) = 0;
                col(j * 2 + 1) = 0;
            }

            float sum = col.sum();
            sum += 1e-30;
            
            renorm[fw][m + sidestep] = srcrenorm + log(sum);
            col *= expf(srcrenorm - renorm[fw][m + sidestep]);
            if (fw && fullwo)
            {
                fwbw[2].col(m + sidestep) *= expf(srcrenorm - renorm[fw][m + sidestep]);
            }
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

float& haplotype::getanyprior(int m) const
{
    return anypriors[m][getindex()];
}

float& haplotype::getnewanyprior(int m) const
{
    return newanypriors[m][getindex()];
}

int basehaps;

template<typename T>
void sortibd2old(T&& probs, int* indices)
{
    double sum = 0;
    for (int x = 0; x < haplotypes.size() / ibdgroupsize; x++)
    {
        indices[x] = x;
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            sum += probs[x * ibdgroupsize * 2 + k];
        }
    }

    std::sort(indices, &indices[haplotypes.size() / ibdgroupsize], [&probs] (int a, int b)
    {
        float asum = 0;
        float bsum = 0;

        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            asum += probs[a * ibdgroupsize * 2 + k];
            bsum += probs[b * ibdgroupsize * 2 + k];
        }
        return asum < bsum;
    });

    bool rerun = false;
    do
    {
    rerun = false;
    for (int x = haplotypes.size() / ibdgroupsize - 1; x >= 0; x--)
    {
        int y = indices[x];
        float ysum = 0;
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            ysum += probs[y * ibdgroupsize * 2 + k];
        }
        if (ysum > 0 && ysum * 3 > (sum + 1e-10) * 1.00001)
        {
            //double diff = sum - ysum * 2;
            double diff = sum - ysum * 3;
            double factor = std::max((diff * 0.5 + ysum) / ysum, 1e-10);
            ///double diff = sum - ysum;
            ///double factor = std::max(diff / ysum, 1e-10);
            for (int k = 0; k < ibdgroupsize * 2; k++)
            {
                probs[y * ibdgroupsize * 2 + k] *= factor;
            }
            if (nodropysum)
            {
               sum -= ysum * (1 - factor);
               if (resetsort) rerun = true;
            }
            //ysum = probs[y * 2] + probs[y * 2 + 1];            
        }
        else if (nodropysum) break;

        if (!nodropysum) sum -= ysum;    
    }
    } while (rerun);
}

template<bool antigroup = ::antigroup, bool alreadyfactor = ibdfactors, size_t N = 1, class T>
void sortibd2b(ArrayXf& probs, std::array<ArrayXXf*, N> first, std::array<ArrayXXf*, N> second, int m, int* indices, ArrayXf& ysums, T&& factors)
{
    double minlevelibd = ::minlevelibd * (mulminibd ? N : 1);
    double sum = 0;
    float max = 0;
    int groupcount = (haplotypes.size() - basehaps) / ibdgroupsize;
    ysums.resize(groupcount);
    factors.resize(groupcount);
    probs.resize(first[0]->col(m).size());
    for (int x = 0; x < basehaps * 2; x++)
    {
        probs[x] = first[0]->col(m)[x] * second[0]->col(m)[x];
        sum += probs[x];
    }
    double basesum = sum;
    for (int x = 0; x < groupcount; x++)
    {
        if (!alreadyfactor) indices[x] = 0;
        double ysum = 0;
        double firstsum = 0;
        double secondsum = 0;
        if (!alreadyfactor) factors[x] = 1.0f;
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            float firstprob = first[0]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k];
            float secondprob = second[0]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k];
            for (int i = 1; i < N; i++)
            {
                if (!ibdmax)
                {
                    firstprob += first[i]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k];
                    secondprob += second[i]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k];
                }
                else
                {
                    firstprob = std::max(firstprob, first[i]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k]);
                    secondprob = std::max(secondprob, second[i]->col(m)[basehaps * 2 + x * ibdgroupsize * 2 + k]);
                }
            }
            float prob = firstprob * secondprob;
            probs[basehaps * 2 + x * ibdgroupsize * 2 + k] = prob * (alreadyfactor ? factors[x] : 1);
            if (!alreadyfactor)
            {
                if (antigroup)
                {
                    if (!ibdmax2)
                    {
                        ysum += prob;
                    }
                    else
                    {
                        ysum = std::max<double>(ysum, prob);
                    }
                }
                else
                {
                    if (!groupmax)
                    {
                        firstsum += firstprob;
                        secondsum += secondprob;
                    }
                    else
                    {
                        firstsum = std::max<double>(firstsum, firstprob);
                        secondsum = std::max<double>(secondsum, secondprob);
                    }
                    /*firstsum += first[x * ibdgroupsize * 2 + k];
                    secondsum += second[x * ibdgroupsize * 2 + k];*/ // TODO
                }
            }
        }
        if (!alreadyfactor)
        {
            if (!antigroup) ysum = firstsum * secondsum;
            sum += ysum;
            ysums[x] = ysum;
            if (adjibd) if (ysum > max) max = ysum;
        }
    }

    if (alreadyfactor) return;

    double orgsum = sum;

    double adjlevelibd = levelibd;
    int indexcount = !adjibd;

    if (adjibd)
    {
        adjlevelibd = std::max(adjlevelibd, declstep * sum / max + minlevelibd - declstep);
        do
        {
            indexcount++;
            double prelevel = adjlevelibd;
            double lim = sum / adjlevelibd;
            int count = 0;
            for (int x = 0; x < groupcount; x++)
            {
                float ysum = ysums[x];

                if (ysum >= lim)
                {
                    count++;
                    if (readj && tightreadj)
                    {
                        double newlevel = declstep * count + minlevelibd - declstep;
                        if (newlevel > adjlevelibd)
                        {
                            adjlevelibd = newlevel;
                            lim = sum / adjlevelibd;
                        }
                    }
                    if (preindex) indices[x] = indexcount;
                }
            }

            adjlevelibd = std::max<double>(minlevelibd, declstep * count + minlevelibd - declstep);
            if (readj && adjlevelibd == prelevel) break;
        } while (readj);        
    }


    bool rerun = false;
    double ibdlevel = adjlevelibd;
    int last = -1;
    do
    {
    rerun = false;
    double sum2 = 0;
    bool anyunchanged = false;
    double fixed = 0;
    for (int x = 0; x < groupcount; x++)
    {
        //if (x == last) continue;
        int y = x;
        float ysum = ysums[x] * factors[x];
        
        if (ibdlevel >= groupcount / 2 - 1)
        {
            ibdlevel = groupcount / 2 - 1;
        }
        else
        {
            anyunchanged = true;
        }

        double nowibdlevel = indices[x] == 0|| ibddecl ? ibdlevel : indices[x];
        double floorlevel = (max * ibdfloor + (orgsumw ? orgsum : 1) * ibdeps);
        if (ysum > floorlevel * 1.00001 && ysum * nowibdlevel > (sum + (orgsumw ? orgsum : 1) * ibdeps) * 1.00001)
        {
            last = x;
            //double diff = sum - ysum * 2;
            double diff = sum - ysum;
            double factor = std::max((diff / (nowibdlevel - 1)), floorlevel) / ysum;

            if (progribd && indices[x] < indexcount)
            {
                double oldibdlevel = ibdlevel;
                ibdlevel += declstep;
                if (dofixed)
                {
                    if (ibdlevel >= groupcount / 2 - 1)
                    {
                        ibdlevel = groupcount / 2 - 1;
                    }
                    if (ibddecl && oldibdlevel != ibdlevel)
                    {
                        double oldfixed = fixed;
                        fixed *= oldibdlevel / ibdlevel;
                        sum -= oldfixed - fixed;
                        sum2 -= oldfixed - fixed;
                    }
                }
            }
            indices[x] = indexcount; // breaking !ibddecl
            ///double diff = sum - ysum;
            ///double factor = std::max(diff / ysum, 1e-10);
            factors[x] *= factor;
            if (nodropysum)
            {
               sum -= ysum * (1 - factor);
               if (resetsort) rerun = true;
            }
            //ysum = probs[y * 2] + probs[y * 2 + 1];
        }       
        if (nodropysum) sum2 += ysums[x] * factors[x];
        if (!nodropysum) sum -= ysum;    
    }
    if (nodropysum) sum = sum2 + basesum;
    rerun &= anyunchanged;
    } while (rerun);

    for (int x = 0; x < groupcount; x++)
    {
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            probs[basehaps * 2 + x * ibdgroupsize * 2 + k] *= factors[x];
        }
    }
    
}


template<typename T>
void sortibd2(T&& probs, int* indices)
{
    double sum = 0;
    float max = 0;
    for (int x = 0; x < haplotypes.size() / ibdgroupsize; x++)
    {
        indices[x] = 0;
        float ysum = 0;
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            if (adjibd) ysum += probs[x * ibdgroupsize * 2 + k];
            sum += probs[x * ibdgroupsize * 2 + k];
        }
        if (adjibd) if (ysum > max) max = ysum;
    }

    double orgsum = sum;

    double adjlevelibd = levelibd;
    int indexcount = !adjibd;

    if (adjibd)
    {
        adjlevelibd = std::max(adjlevelibd, declstep * sum / max + minlevelibd - declstep);
        do
        {
            indexcount++;
            double prelevel = adjlevelibd;
            double lim = sum / adjlevelibd;
            int count = 0;
            for (int x = 0; x < haplotypes.size() / ibdgroupsize; x++)
            {
                float ysum = 0;
                for (int k = 0; k < ibdgroupsize * 2; k++)
                {
                    ysum += probs[x * ibdgroupsize * 2 + k];
                }         
                if (ysum >= lim)
                {
                    count++;
                    if (readj && tightreadj)
                    {
                        double newlevel = declstep * count + minlevelibd - declstep;
                        if (newlevel > adjlevelibd)
                        {
                            adjlevelibd = newlevel;
                            lim = sum / adjlevelibd;
                        }
                    }
                    if (preindex) indices[x] = indexcount;
                }
            }

            adjlevelibd = std::max<double>(minlevelibd, declstep * count + minlevelibd - declstep);
            if (readj && adjlevelibd == prelevel) break;
        } while (readj);        
    }


    bool rerun = false;
    double ibdlevel = adjlevelibd;
    int last = -1;
    do
    {
    rerun = false;
    double sum2 = 0;
    bool anyunchanged = false;
    double fixed = 0;
    for (int x = 0; x < haplotypes.size() / ibdgroupsize; x++)
    {
        //if (x == last) continue;
        int y = x;
        float ysum = 0;
        for (int k = 0; k < ibdgroupsize * 2; k++)
        {
            ysum += probs[y * ibdgroupsize * 2 + k];
            if (nodropysum) sum2 += probs[y * ibdgroupsize * 2 + k];
        }
        if (ibdlevel >= haplotypes.size() / ibdgroupsize / 2 - 1)
        {
            ibdlevel = haplotypes.size() / ibdgroupsize / 2 - 1;
        }
        else
        {
            anyunchanged = true;
        }

        double nowibdlevel = indices[x] == 0|| ibddecl ? ibdlevel : indices[x];
        double floorlevel = (max * ibdfloor + (orgsumw ? orgsum : 1) * ibdeps);
        if (ysum > floorlevel * 1.00001 && ysum * nowibdlevel > (sum + (orgsumw ? orgsum : 1) * ibdeps) * 1.00001)
        {
            last = x;
            //double diff = sum - ysum * 2;
            double diff = sum - ysum;
            double factor = std::max((diff / (nowibdlevel - 1)), floorlevel) / ysum;

            if (progribd && indices[x] < indexcount)
            {
                double oldibdlevel = ibdlevel;
                ibdlevel += declstep;
                if (dofixed)
                {
                    if (ibdlevel >= haplotypes.size() / ibdgroupsize / 2 - 1)
                    {
                        ibdlevel = haplotypes.size() / ibdgroupsize / 2 - 1;
                    }
                    if (ibddecl && oldibdlevel != ibdlevel)
                    {
                        double oldfixed = fixed;
                        fixed *= oldibdlevel / ibdlevel;
                        sum -= oldfixed - fixed;
                        sum2 -= oldfixed - fixed;
                    }
                }
            }
            indices[x] = indexcount; // breaking !ibddecl
            ///double diff = sum - ysum;
            ///double factor = std::max(diff / ysum, 1e-10);
            for (int k = 0; k < ibdgroupsize * 2; k++)
            {
                probs[y * ibdgroupsize * 2 + k] *= factor;
            }
            if (nodropysum)
            {
               sum -= ysum * (1 - factor);
               sum2 -= ysum * (1 - factor);
               if (dofixed && ibddecl && progribd) fixed += ysum * factor;
               if (resetsort) rerun = true;
            }
            //ysum = probs[y * 2] + probs[y * 2 + 1];            
        }

        if (!nodropysum) sum -= ysum;    
    }
    if (nodropysum) sum = sum2;
    rerun &= anyunchanged;
    } while (rerun);
}


template<class column> void doemit(column& c, float& anyprior, genprob& prior, int marker, int* indices)
{
    vector<genprob>& ourPrior = priors[marker];
    vector<float>& ourAnyPrior = anypriors[marker];
    #pragma ivdep
    for (int i = 0; i < haplotypes.size(); i++)
    {
        float old = c[i * 2] + c[i * 2 + 1];
        float sum = 0;
        for (int j = 0; j < 2; j++)
        {
            float val = (anyprior ? prior[j] : 1.0f) * ourPrior[i][j];
            sum += val /** val*/;
        }
        /*sum = 1;
        float norm = 0;
        for (int j = 0; j < 2; j++)
        {
            norm += ourPrior[i][j] * ourPrior[i][j];
        }
        sum *= norm;*/        
        sum = 1;
        for (int j = 0; j < 2; j++)
        {
            float val = sum;
            val *= (anyprior ? prior[j] : 1.0f) * ourPrior[i][j];
            
            float anyPriorW = /*anyprior * */ourAnyPrior[i] ? 1.0f : 0.0f;
            val *= anyPriorW;
            val += 0.5f * (1.0f - anyPriorW) * (anyprior ? prior[j] : 1.0f);
            c[i * 2 + j] = old * val * (scaleanyprior ? anyprior : 1.0f) * (scaleanypriorw ? ourAnyPrior[i] : 1.0f);
//        if (val < 0 || val > 1) printf("%f\n", val);
        }
    }

    if (sortemit)
    {
        sortibd2(c, indices);
    }
}

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d, int index)
{
    // Careful! c and c2 might coincide
    float dist = (themap.chromposes[marker + d] - themap.chromposes[marker]) * d * -0.02 * Ne;
    float nonrec = expf(dist);
    int actualSize = haplotypes.size();
    if (onlyref) actualSize = basehaps;
    if (halfpar) actualSize -= basehaps / 2;
    float rec = -expm1f(dist) / actualSize;
    float sum = c.sum();
    float subsum = 0;
    float subunc = 0;
    float certf = 0;
    float partsum = 0;
    int prevbase = -1;
    for (int i = 0; i < haplotypes.size(); i++)
    {
        if (killibd2trans)
        {
            int base = i / ploidy * ploidy;
            if (base != prevbase)
            {
                subsum = c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1)).sum();
                if (killibd2unc)
                {                    
                    partsum = (c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1, 2)) + c(Eigen::seq(base * 2 + 1, (base + ploidy) * 2 - 1, 2))).square().sum();
                    subunc = partsum / (subsum * subsum + 1e-30f);
                }
                if (genounc)
                {
                    float term1 = c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1, 2)).sum();
                    float term2 = c(Eigen::seq(base * 2 + 1, (base + ploidy) * 2 - 1, 2)).sum();
                    certf = (term1 * term1 + term2 * term2) / (subsum * subsum + 1e-30f);
                }
                if (sepgenounc)
                {
                    certf = c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1, 2)).square().sum() + c(Eigen::seq(base * 2 + 1, (base + ploidy) * 2 - 1, 2)).square().sum();
                    certf /= partsum + 1e-30f;
                }

                prevbase = base;
            }
        }
        float old = c[i * 2] + c[i * 2 + 1];
        bool filter = (onlyref && i >= basehaps);
        if (halfpar) filter |= (i < basehaps) && (((index - basehaps) % ploidy < ploidy / 2) ^ (i < basehaps / 2));
        if (filterrefs && i < basehaps)
        {
            bool ok = false;
            if (haplotypes[index].allowedrefs)
            {
                for (int ref : *haplotypes[index].allowedrefs)
                {
                    if (ref == i)
                    {
                        ok = true;
                        break;
                    }
                }
            }
            filter |= !ok;
        }
        for (int j = 0; j < 2; j++)
        {
            c2[i * 2 + j] = (1.0f - (filter ? filterlevel : 0)) * ((old * (certf + (1 - certf) * (killibd2trans ? (killibd2unc ? subunc : std::max(old, subsum - old) / (subsum + 1e-30f)) : 1 ))) * nonrec + sum * rec);
        }
    }
}

std::mt19937 rng;

struct individ
{
    vector<int> genotypes;
    vector<array<int, 2>> reads;
    array<int, 8> allowedrefs;

    individ()
    {
        for (auto& allowed : allowedrefs)
        {
            allowed = -1;
        }
    }
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
    std::uniform_real_distribution<float> distribution(-offsmagn, offsmagn);

    for (int j = 0; j < ploidy; j++)
    {        
        haplotypes[index + j].posterior.resize(genotypes.size());
        haplotypes[index + j].posteriorwo.resize(genotypes.size());
        haplotypes[index + j].offset.resize(genotypes.size());
        haplotypes[index + j].sim.resize(genotypes.size());
        haplotypes[index + j].crosssim.resize(genotypes.size());
        haplotypes[index + j].desired.resize(genotypes.size());
        haplotypes[index + j].momentum.resize(genotypes.size());
        if (allowedrefs[0] != -1) haplotypes[index + j].allowedrefs = &allowedrefs;
        bool first = true;
        int firstatall = -1;
        for (int i = 0; i < genotypes.size(); i++)
        {
            //if (ourmap.otherepses[i] > 0.10) genotypes[i] = -1;
            double genotype = genotypes[i];
            bool any = false;
            if (genotype < 0)
            {
                int readsum = reads[i][0] + reads[i][1];
                if (readsum || updallpriors)
                {
                    genotype = (reads[i][1] + (ploidy - 1) * 0.5) / (readsum + ploidy - 1) * ploidy;
                    haplotypes[index + j].getanyprior(i) = std::max(0.5f, 1.0f - powf(powf(0.5, 1.0f / ploidy), readsum));
                }
                else
                    haplotypes[index + j].getanyprior(i) = false;
                any = reads[i][0] && reads[i][1];
            }
            else
            {
                haplotypes[index + j].getanyprior(i) = true;
                any = genotype >= 1 && genotype <= ploidy - 1;
            }

            haplotypes[index + j].offset[i] = distribution(rng);
            haplotypes[index + j].momentum[i] = 0;

            if (genotype >= 0)
            {
                if (firstatall == -1) firstatall = i;
                
                float val = std::clamp<float>((genotype / 1.0f / ploidy) * (1.0f - haplotypes[index + j].offset[i]), updeps, 1 - updeps);
                haplotypes[index + j].getprior(i)[0] = 1.0f - val;
                haplotypes[index + j].getprior(i)[1] = val;                
                if (fixatone)
                {                
                if (!nozeroone && first && any && (j == 0 || j == ploidy - 1))
                {
                    first = false;

                    float a = 1.0f;
                    float b = 0.f;
                    if (j == ploidy - 1)
                    {
                        std::swap(a, b);
                    }

                    haplotypes[index + j].getprior(i)[0] = a;
                    haplotypes[index + j].getprior(i)[1] = b;
                }
                }
            }
            else
            {
                haplotypes[index + j].getprior(i)[0] = 0;
                haplotypes[index + j].getprior(i)[1] = 0;    
            }
        }

        if (first && !j && firstatall >= 0)
        {
            bool val = haplotypes[index + j].getprior(firstatall)[0] < 0.5;

            if (fixatone)
            {
            haplotypes[index + j].getprior(firstatall)[0] = !val;
            haplotypes[index + j].getprior(firstatall)[1] = val;
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
// TODO LESS MEMORY
    //double scores[haplotypes[index].fwbw[0].cols()][permcount];    
    vector<array<double, permcount> > scores;
    scores.resize(haplotypes[index].fwbw[0].cols());
    int indices[haplotypes.size() / ibdgroupsize];
    ArrayXf probs[ploidy], ysums, factors;
//    #pragma omp parallel for schedule(guided, 100), private(indices, probs, ysums, factors), shared(scores), num_threads(ploidy * 2)
    #pragma omp taskloop num_tasks(ploidy * 2), private(indices, probs, ysums, factors), shared(scores)
    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        bool first = true;
        double firstthisscore = 0;
        double firstagnscore = 0;
        int firstpow2 = 0;
        float sims[ploidy][ploidy];
        float corrs[ploidy];

        if (ibdfactors)
        {
            for (int k = 0; k < (antifactors ? ploidy : 1); k++)
            {
                constexpr int count = antifactors ? 1 : ploidy;
                std::array<ArrayXXf*, count> first;
                std::array<ArrayXXf*, count> firstwo;
                std::array<ArrayXXf*, count> second;

                for (int j = 0; j < count; j++)        
                {                    
                    first[j] = &haplotypes[index + k + j].fwbw[1 + powofactors];
                    firstwo[j] = &haplotypes[index + k + j].fwbw[1 + fullwo];
                    second[j] = &haplotypes[index + k + j].fwbw[0];
                }
                sortibd2b<antigroup, false, count>(probs[0], first, second, m, indices, ysums, haplotypes[index + k].fwbw[2 + fullwo].col(m));
                if (nonsimfactor) sortibd2b<antigroup, false, count>(probs[0], firstwo, second, m, indices, ysums, haplotypes[index + k].fwbw[2 + fullwo + nonsimfactor].col(m));
            }
        }

        #pragma ivdep
        for (int j = 0; j < ploidy; j++)        
        {                    
            if (!expklsim)
            {
                if (sortsim && groupibd2)
                {
                    sortibd2b(probs[j], {&haplotypes[index + j].fwbw[1 + simpowo]}, {&haplotypes[index + j].fwbw[0]}, m, indices, ysums, haplotypes[index + (antifactors ? j : 0)].fwbw[2 + fullwo + nonsimfactor].col(m));
                }
                else
                {
                    probs[j] = haplotypes[index + j].fwbw[1 + simpowo].col(m) * haplotypes[index + j].fwbw[0].col(m);
                    if (sortsim)
                    {
                        sortibd2(probs[j], indices);
                    }
                }
                corrs[j] = 1.0f / (sqrt((probs[j] * probs[j]).sum()) + 1e-30f);
            }
            else
            {
                probs[j] = haplotypes[index + j].fwbw[1 + simpowo].col(m) * haplotypes[index + j].fwbw[0].col(m) + 1e-30f;
                if (sortsim) sortibd2(probs[j], indices);
                corrs[j] = 1.0f / probs[j].sum();
            }
        }

        for (int j = 0; j < ploidy; j++)        
        {
            float& sim = haplotypes[index + j].sim[m];
            sim = -1000.f;

            for (int k = 0; k < ploidy; k++)
            {                
                if (k == j)
                {
                    continue;
                }
                if (!expklsim)
                {
                    haplotypes[index + j].crosssim[m][k] = std::clamp<float>((probs[j] * probs[k]).sum() * corrs[j] * corrs[k], 0, 1);
                }
                else
                {
                    haplotypes[index + j].crosssim[m][k] = exp(-(probs[j] * (probs[j].log() - probs[k].log())).sum() * corrs[j] - log(corrs[j]) + log(corrs[k]));
                }
                sim = std::max(sim, haplotypes[index + j].crosssim[m][k]); 
            }
        }

        int pow2s[ploidy][ploidy];
        double singlescores[ploidy][ploidy];

        for (int j = 0; j < ploidy; j++)
        {
            for (int k = 0; k < ploidy; k++)
            {
                double sumterm = 0;
                if (sortflip && groupibd2)
                {
                    if (!vetoflip && !crossibd)
                    {
                        sortibd2b<antigflip>(probs[k], {&haplotypes[index + j].fwbw[1]}, {&haplotypes[index + k].fwbw[0]}, m, indices, ysums, haplotypes[index + (antifactors ? j : 0)].fwbw[2 + fullwo].col(m));
                    }
                    else
                    {
                        sortibd2b<antigflip, false>(probs[k], {&haplotypes[index + j].fwbw[1]}, {&haplotypes[index + k].fwbw[0]}, m, indices, ysums, factors);
                    }
                }
                else
                {
                    probs[k] = haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + k].fwbw[0].col(m);
                    if (sortflip) sortibd2(probs[k], indices);
                }
                // TODO returnera summan så vi har den
                for (int i = 0; i < haplotypes.size() * 2; i++)
                {
                    /*double a = haplotypes[index + j].fwbw[1].col(m)(i);
                    double b = haplotypes[index + k].fwbw[0].col(m)(i);*/

                    sumterm += probs[k](i);
                }
                if (vetoflip && j != k)
                {
                    for (int z : {j, k})
                    {
                        double sumterm2  = 0;
                        sortibd2b<antigflip>(probs[k], {&haplotypes[index + j].fwbw[1]}, {&haplotypes[index + k].fwbw[0]}, m, indices, ysums, haplotypes[index + (antifactors ? z : 0)].fwbw[2 + fullwo].col(m));
                        for (int i = 0; i < haplotypes.size() * 2; i++)
                        {
                            /*double a = haplotypes[index + j].fwbw[1].col(m)(i);
                            double b = haplotypes[index + k].fwbw[0].col(m)(i);*/

                            sumterm2 += probs[k](i);
                        }
                        if (sumterm2 < sumterm) sumterm = sumterm2;
                    }
                }
                singlescores[j][k] = log(frexp(sumterm, &pow2s[j][k]));
            }
        }

        for (int p = permcount - 1; p >= 0; p--)
        {
            array<int, ploidy> perm;
            bool badperm = !getploidyperm(p, perm);
            if (!badperm && (halfpar || true))
            {
                for (int i = 0; i < ploidy; i++)
                {
                    if ((perm[i] < ploidy / 2) ^ (i < ploidy / 2)) badperm = true;
                }
            }
            if (badperm)
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
            double sumagn = 0;
            int pow2 = 0;
            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                //sum += haplotypes[index + j].renorm[1][m];
                auto sumuptoploidy = [this](const auto& vector1, const auto& vector2, int baseindex)
                {
                    return (vector1.reshaped(haplotypes.size() / ploidy / 2, ploidy * 2)(Eigen::all, Eigen::seq(baseindex, 2 * ploidy - 1, 2)).rowwise().sum().eval() * 
                    vector2.reshaped(haplotypes.size() / ploidy / 2, ploidy * 2)(Eigen::all, Eigen::seq(baseindex, 2 * ploidy - 1, 2)).rowwise().sum().eval()).sum();
                };
                if (agnosticflip)
                {
                    float sumagnterm = 0;
                    for (int i = 0; i < haplotypes.size() / ploidy; i++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            float terms[2] = {0};
                            for (int n = 0; n < ploidy; n++)
                            {
                                for (int z = 0; z < 2; z++)
                                { 
                                    int subindex = ((i * ploidy + n) * 2 + k);
                                    int nowindex = k ? index + j : index + perm[j];
                                    terms[z] += haplotypes[nowindex].fwbw[k].col(m)(subindex);
                                }
                            }
                            sumagnterm += terms[0] * terms[1];
                        }
                    }
                    /*sumuptoploidy(haplotypes[index + j].fwbw[1].col(m), haplotypes[index + perm[j]].fwbw[0].col(m), 0);
                    sumagnterm += sumuptoploidy(haplotypes[index + j].fwbw[1].col(m), haplotypes[index + perm[j]].fwbw[0].col(m), 1);*/
                    sumagn += log(sumagnterm + 1e-30);
                }
                if (!agnosticflip || bothagn)
                {
                    //sum += log((haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + perm[j]].fwbw[0].col(m)).sum() + 1e-30);
                    double sumterm = 0;
                    /*for (double val : (haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + perm[j]].fwbw[0].col(m)))
                    {
                        sumterm += val;
                    }*/

                    //sum += log(sumterm);
                    sum += singlescores[j][perm[j]];
                    pow2 += pow2s[j][perm[j]];
                }
                //sum += haplotypes[index + perm[j]].renorm[0][m];
            }

            if (agnosticflip && !bothagn)
            {
                sum = sumagn;
            }

            if (first)
            {
                if (m == 0)
                {
                    double firstscore = sum;
                    for (int j = 0; j < ploidy; j++)
                    {
                        firstscore += haplotypes[index + j].renorm[1][m];   
                        firstscore += haplotypes[index + j].renorm[0][m];
                    }
                    firstscore += log(2) * pow2;
                    #pragma omp atomic
                    likelihood += firstscore;
                }
                sum += 0.001;
                firstpow2 = pow2;
                firstthisscore = sum;
                firstagnscore = sumagn;
                first = false;
            }

            //if (index == 16) printf("Flip: %d %d %d %f\n", index, m, p, sum);
            sum += log(2) * (pow2 - firstpow2);
            scores[m][p] = sum - firstthisscore;
            if (bothagn && sumagn < firstagnscore) scores[m][p] = sumagn - firstagnscore;
        }
    }

    #pragma omp taskwait

    int bestmarker = 0;
    int bestp = 0;
    double bestscore = -1.1e30f;
    double realbestscore = -1.1e30f;    
    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        for (int p = permcount - 1; p >= 0; p--)
        {
            double sum = scores[m][p] * (noiseflippos ? std::uniform_real_distribution<double>(0, 1)(rng) : 1);
            if (sum < -1e30f) continue;

            if (!randflippos || bestscore < 0)
            {
                if (sum > bestscore)
                {
                    bestscore = sum;
                    realbestscore = sum;
                    bestp = p;
                    bestmarker = m;
                }
            }
            else if (sum >= 0)
            {
                if (!std::bernoulli_distribution(1 / (exp(sum - bestscore) + 1))(rng))
                {                
                    bestp = p;
                    bestmarker = m;                    
                }
                if (sum > realbestscore)
                {
                    realbestscore = sum;
                }
                bestscore += log(exp(sum - bestscore) + 1);
            }
        }
    }    

    return {bestmarker, bestp, realbestscore};
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

    straight |= std::bernoulli_distribution(rngflwght ? flipscale / (exp(bestscore - flipdrag) + 1) : 0.5)(rng) || neverflip;
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
        array<genprob, ploidy> posteriorwo;
        array<float, ploidy> offset;
        array<float, ploidy> sim;
        array<float, ploidy> desired;
        array<float, ploidy> momentum;
        array<float, ploidy> anyprior;
        array<float, ploidy> newanyprior;
        for (int i = 0; i < haplotypes[index].posterior.size(); i++)
        {
            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                prior[j] = haplotypes[index + j].getprior(i);
                posterior[j] = haplotypes[index + j].posterior[i];
                posteriorwo[j] = haplotypes[index + j].posteriorwo[i];
                offset[j] = haplotypes[index + j].offset[i];
                sim[j] = haplotypes[index + j].sim[i];
                desired[j] = haplotypes[index + j].desired[i];
                momentum[j] = haplotypes[index + j].momentum[i] * (1.0f - killflipm);
                if (dovar || wounc)
                {
                    anyprior[j] = haplotypes[index + j].getanyprior(i);
                    newanyprior[j] = haplotypes[index + j].getnewanyprior(i);
                }
            }

            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                int permval = i > bestmarker ? perm[j] : j;
                haplotypes[index + j].getnewprior(i) = prior[permval];
                haplotypes[index + j].posterior[i] = posterior[permval];
                haplotypes[index + j].posteriorwo[i] = posteriorwo[permval];
                haplotypes[index + j].offset[i] = offset[permval];
                haplotypes[index + j].sim[i] = sim[permval];
                haplotypes[index + j].desired[i] = desired[permval];
                haplotypes[index + j].momentum[i] = momentum[permval];                
                if (dovar || wounc)
                {
                    haplotypes[index + j].getanyprior(i) = anyprior[permval];
                    haplotypes[index + j].getnewanyprior(i) = newanyprior[permval];
                }
            }
        }
    }

    return !straight;
}

void individ::doposteriorhaplotypes(int index)
{
    ArrayXf probs, probswo, weights, ysums;
    int indices[haplotypes.size() / ibdgroupsize];
    #pragma omp taskloop num_tasks(ploidy * 2), private(probs, probswo, weights, indices, ysums), collapse(2)
//#pragma omp parallel for private(probs, probswo, weights, indices, ysums), collapse(2), schedule(dynamic, 100), num_threads(ploidy * 2)
    for (int j = 0; j < ploidy; j++)
    {
        for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
        {
            if (ibd2sort && groupibd2)
            {
                sortibd2b(probs, {&haplotypes[index + j].fwbw[1]}, {&haplotypes[index + j].fwbw[0]}, m, indices, ysums, haplotypes[index + (ibdfactors && !antifactors ? 0 : j)].fwbw[2 + fullwo].col(m));
            }
            else
            {
                probs = haplotypes[index + j].fwbw[1].col(m) * haplotypes[index + j].fwbw[0].col(m);
                if (ibd2sort)
                {
                    sortibd2(probs, indices);
                }
            }

            if (fullwo)
            {
                if (ibd2sort && groupibd2)
                {
                    sortibd2b(probswo, {&haplotypes[index + j].fwbw[2]}, {&haplotypes[index + j].fwbw[0]}, m, indices, ysums, haplotypes[index + (ibdfactors && !antifactors ? 0 : j)].fwbw[2 + fullwo + nonsimfactor].col(m));
                }
                else
                {
                    probswo = haplotypes[index + j].fwbw[2].col(m) * haplotypes[index + j].fwbw[0].col(m);
                    if (ibd2sort)
                    {
                        sortibd2(probswo, indices);
                    }
                }
            }
                        
            //if (!antidiffp)
            double var = 0;
            if (!advpostsum)
            for (int z = 0; z < 2; z++)
            {            
                haplotypes[index + j].posterior[m][z] = probs(Eigen::seq(z, haplotypes.size() * 2 - 1, 2)).sum();
            }
            else
            {
                double vals[std::max(advpostsum * 2, 3)] = {0};
                for (int i = 0; i < probs.size(); i += 2)
                {
                    double part = 0;
                    for (int z = 0; z < 2; z++)
                    {
                        vals[z] += probs(i + z);
                        //if (dovar) part += probs(i + z);
                    }
                    //if (dovar && part) vals[2] += vals[z] * vals[z] / part;
                }
                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posterior[m][z] = vals[z];
                }
                //if (dovar) var += vals[2];
            }
            //if (antidiffp)
            /*if (fullwo)
            {
                for (int z = 0; z < 2; z++)
                    haplotypes[index + j].posteriorwo[m][z] = (haplotypes[index + j].fwbw[2].col(m) * haplotypes[index + j].fwbw[0].col(m))(Eigen::seq(z, haplotypes.size() * 2 - 1, 2)).sum();
            }
            else*/
            if (selfposteriorwo)
            {
                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posteriorwo[m][z] = 0;
                }
                auto& woprobs = fullwo ? probswo : probs;
                double wovals[std::max(advpostsum * 2, 3)] = {0};
                for (int i = 0; i < woprobs.size(); i += 2)
                {
                    /*float weight = priors[m][i / 2][0] * priors[m][i / 2][0] + priors[m][i / 2][1] * priors[m][i / 2][1];
                    float weight2 = probs(i) * priors[m][i / 2][0] + probs(i + 1) * priors[m][i / 2][1];
                    float sum2 = probs(i) * probs(i) + probs(i + 1) * probs(i + 1);

                    weight = weight2 * weight2 / (weight +1e-30f) / (sum2 + 1e-30f);*/                    
                    float sum2 = woprobs(i) + woprobs(i + 1);
                    float prior2 = priors[m][i / 2][0] + priors[m][i / 2][1];                
                    if (prior2)
                    for (int z = 0; z < 2; z++)
                    {
                        if (advpostsum)
                        {
                            wovals[z] += sum2 * priors[m][i / 2][z];
                        }
                        else
                        {
                            haplotypes[index + j].posteriorwo[m][z] += sum2 * priors[m][i / 2][z];
                        }
                    }
                    else if (fillinmissingwo)
                    for (int z = 0; z < 2; z++)
                    {
                        assert(false);
                        haplotypes[index + j].posteriorwo[m][z] += sum2 * (neutralmissing ? 0.5f : priors[m][index + j][z]);
                    }
                }

                if (advpostsum)
                for (int z = 0; z < 2; z++)
                {                    
                    haplotypes[index + j].posteriorwo[m][z] = wovals[z];
                }
            }
            if (false)
            for (int i = 0; i < probs.size(); i += 2)
            {
                /*float weight = priors[m][i / 2][0] * priors[m][i / 2][0] + priors[m][i / 2][1] * priors[m][i / 2][1];
                float weight2 = probs(i) * priors[m][i / 2][0] + probs(i + 1) * priors[m][i / 2][1];
                float sum2 = probs(i) * probs(i) + probs(i + 1) * probs(i + 1);

                weight = weight2 * weight2 / (weight +1e-30f) / (sum2 + 1e-30f);*/
                float weight = probs(i) * probs(i) + probs(i + 1) * probs(i + 1);
                float sum2 = probs(i) + probs(i + 1);
                weight /= sum2 * sum2 + 1e-30f;
                bool skip = false;
                if (!weight) skip = true;
                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posterior[m][z] += probs(i + z) * (!skip ? weight : 0.5f);
                    //haplotypes[index + j].posterior[m][z] += probs(i + z) * (!skip ? priors[m][i / 2][z] * priors[m][i / 2][z] / weight : 0.5f);
                }
            }


            float sum = 1e-30f;
            float sumwo = 1e-30f;
            for (int z = 0; z < 2; z++)
            {
                sum += haplotypes[index + j].posterior[m][z];
                if (selfposteriorwo) sumwo += haplotypes[index + j].posteriorwo[m][z];
            }
            if (sum < 1e-12f) printf("HEJ %d %d %g\n", index + j, m, sum);
            if (selfposteriorwo && sumwo < 1e-12f) printf("HEJWO %d %d %g\n", index + j, m, sumwo);
            double origsum = sum;
            sum = 1 / sum;
            if (selfposteriorwo) sumwo = 1 / sumwo;
            auto& priors = haplotypes[index + j].getprior(m);
            double uncvar = 0;
            for (int z = 0; z < 2; z++)
            {
                haplotypes[index + j].posterior[m][z] *= sum;
                if (dovar)
                {
                    if (!dovarunc || scaleunc)
                    {
                        var += haplotypes[index + j].posterior[m][z] * haplotypes[index + j].posterior[m][z];                        
                    }
                    if (!onlyvar)
                    {
                        uncvar += haplotypes[index + j].posterior[m][z] * priors[z];
                    }
                }
                if (selfposteriorwo) haplotypes[index + j].posteriorwo[m][z] *= sumwo;
            }
            if (haplotypes[index + j].getanyprior(m))
            {
            if (dovar)
            {
                if (scaleunc) uncvar = std::max((((uncvar / var) - uncshift) / (1 - uncshift)), 1e-30);
                    haplotypes[index + j].getnewanyprior(m) = std::max(1e-30, dovarunc ? uncvar : (onlyvar ? var : (1 - var + uncvar)));
            }
            if (wounc)
            {
                haplotypes[index + j].getnewanyprior(m) = origsum * sumwo;
                }
                if ((dovar || wounc) && uncnogeno && (genotypes[m] == -1 && !reads[m][0] && !reads[m][1]))
                {
                    haplotypes[index + j].getnewanyprior(m) *= 0.5;
                }
            }
        }
    }
}

void individ::nudgehaplotypes(int index)
{
    auto updatenewpriors = [this, index] (int i, array<ratiotype, ploidy>& ratio)
    {
        array<double, ploidy> val, step, midpoints, origmidpoints;
        double abssum = 0;
        double plainsum = 0;
        double maxcentered = 0.0;
        bool poshz = (!reads[i][0] || !reads[i][1]) && (genotypes[i] == -1 || genotypes[i] % ploidy == 0);

        auto doextremis = [this, &ratio, index, i]
        {
            bool extremes[ploidy] = {0};
            if (extremsim && !burnin)
            {
                for (int z = 0; z < 2; z++)
                {
                    if (!reads[i][z] && (genotypes[i] == -1 || genotypes[i] == (!z) * ploidy)) continue;
    
                    int extremis = 0;
                    for (int m = 1; m < ploidy; m++)
                    {
                        float diff = ratio[m] - ratio[extremis];
                        if (z) diff *= -1;
                        if (diff > 0) extremis = m;
                    }
                    float sim = haplotypes[index + extremis].sim[i] * dampextreme * (domarkeps ? (1.0f - ourmap.otherepses[i]) : 1.0f);
                    extremes[extremis] = true;
                    ratio[extremis] = ratio[extremis] * (1 - sim) + !z * sim;                    
                }
            }
            if (highsimpost)
            {
                for (int m = 0; m < ploidy; m++)
                {
                    if (extremes[m]) continue;
                    float sim = haplotypes[index + m].sim[i];
                    ratio[m] = ratio[m] * (1 - sim) + haplotypes[index + m].posterior[i][!antihisim] * sim;
                }
            }
        };
        if (preextremis) doextremis();

        for (int m = 0; m < ploidy; m++)
        {
            double midpoint = ratio[m];
            double midpoint1m = 1 - midpoint;
            midpoints[m] = log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));
            origmidpoints[m] = midpoints[m];
            if (priororig)
            {
                midpoint = haplotypes[index + m].getprior(i)[0];
                midpoint1m = haplotypes[index + m].getprior(i)[1];
                if (midpoint11m) midpoint1m = 1 - midpoint;
                origmidpoints[m] = log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));
            }
            if (priorpowo || mixorig)
            {
                midpoint = haplotypes[index + m].posteriorwo[i][0];
                midpoint1m = haplotypes[index + m].posteriorwo[i][1];
                if (midpoint11m) midpoint1m = 1 - midpoint;
                if (midpoint11m) midpoint = 1 - haplotypes[index + m].posteriorwo[i][1];
                origmidpoints[m] = log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));
                if (mixorig)
                {
                    midpoint = haplotypes[index + m].posterior[i][0];
                    midpoint1m = haplotypes[index + m].posterior[i][1];
                    midpoint1m = 1 - midpoint;
                    origmidpoints[m] = (1 - posteriormix) * origmidpoints[m] + posteriormix * log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));                    
                }
            }
            if (postorig)
            {
                midpoint = haplotypes[index + m].posterior[i][0];
                midpoint1m = haplotypes[index + m].posterior[i][1];
                midpoint1m = 1 - midpoint;
                origmidpoints[m] = log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));
            }
        }

        if (crosssimoffset && !burnin && (!nocshz || poshz))
        {
            for (int m = 0; m < ploidy; m++)
            {
                auto& newpriorsm = haplotypes[index + m].getnewprior(i);
                for (int k = 0; k < ploidy; k++)
                {
                    if (m == k) continue;
                    auto& newpriorsk = haplotypes[index + m].getnewprior(i);
                    if (caponpriors)
                    {
                        if (newpriorsm[0] == 0 || newpriorsk[0] == 0 || newpriorsm[1] == 0 || newpriorsk[1] == 0) continue;
                    }
                    else
                    {
                        if (fabs(origmidpoints[m]) > 30 || fabs(origmidpoints[k]) > 30) continue;
                    }

                    // 1 in numerator implies switched order
                    double diff = plaincsdiff ? 1/(1+exp(origmidpoints[k])) - 1/(1+exp(origmidpoints[m])) : (origmidpoints[m] - origmidpoints[k]);
                    if (selflipcs && (midpoints[m] - midpoints[k]) * diff < 0) diff = -diff;
                    if (rewcsdiff) diff = diff / (1.0001 - haplotypes[index + m].crosssim[i][k]);
                    midpoints[m] += (diff + (diff < 0 ? -1 : 1) * csbump) * (fabs(haplotypes[index + m].offset[i]) + fabs(haplotypes[index + k].offset[i]))
                        * (invcs ? 1.0 / (1.0001 - haplotypes[index + m].crosssim[i][k]) - 1.0 + 0.0001 : haplotypes[index + m].crosssim[i][k])
                        * csscale;                    
                }
               midpoints[m] = std::clamp<double>(midpoints[m], -midpointcap, midpointcap);
               ratio[m] = exp(midpoints[m]) / (1 + exp(midpoints[m]));
            }
        }

        if (!preextremis) doextremis();

        for (int m = 0; m < ploidy; m++)
        {
            auto& priors = haplotypes[index + m].getnewprior(i);
            haplotypes[index + m].desired[i] = ratio[m];

            double midpoint = ratio[m];
            midpoint = log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(1 - midpoint, 1e-10, 1.));
            if (!extremetension) midpoint *= tension + (tensionoffset ? haplotypes[index + m].offset[i] : 0.0f);
            midpoints[m] = midpoint;
            //if (burnin) midpoint *= ploidy;            
            double num = std::clamp<double>(priors[0], 1e-10, 1.);
            double denom = std::clamp<double>(priors[1], 1e-10, 1.);
            if (domaxcentered) maxcentered = std::max(maxcentered, 2.0 * std::min(priors[0], priors[1]));

            val[m] = log(num/denom);
            if (extremetension && val[m] * midpoint > 0 && fabs(midpoint) > fabs(val[m])) midpoint *= tension + (tensionoffset ? haplotypes[index + m].offset[i] : 0.0f);
            if (!simplestep)
            {
                if (logitstep)
                {
                    step[m] = midpoint - val[m];
                }
                else
                if (advstep)
                {
                    double a = exp(midpoint);
                    double b = exp(val[m]);
                    //step[m] = (a * (b + 1) - b * (a + 1)) / ((a + 1) * (b + 1));
                    step[m] = (a - b) / ((a + 1) * (b + 1));
                    if (a < b && step[m] > -minstep) step[m] = -minstep;
                    if (a > b && step[m] < minstep) step[m] = minstep;
                }
                else
                {
                    midpoint = exp(midpoint) / (1 + exp(midpoint));
                    step[m] = 1.0 / (exp(val[m]) + 1) + midpoint - 1.0;
                }
            }
            else
            {
                step[m] = midpoint - val[m];
            }
            if (index == 0 && i == 3) printf("\n %d %d %d %lf %lf\n", index, m, i, val[m], step[m]);

            abssum += fabs(step[m]);
            plainsum += step[m];
        }

        plainsum = fabs(plainsum) + 1e-10f;

        if (!tensionoffset && !crosssimoffset)
        for (int m = 0; m < ploidy; m++)
        {
            for (int n = 0; n < ploidy; n++)
            {
                if (m == n) continue;
                if (val[m] - val[n] < -0.01 && midpoints[n] - val[m] < -0.01 && val[n] - midpoints[m] < -0.01)
                {
                    std::swap(val[m], val[n]);
                    std::swap(midpoints[m], midpoints[n]);
                    std::swap(step[m], step[n]);
                    std::swap(haplotypes[index + m].momentum[i], haplotypes[index + n].momentum[i]);
                }
                //if (val[m] < 0 && step[m] < 0 && val[n] > 0 && step[n] > 0 && haplotypes[index + m].offset[i] > haplotypes[index + n].offset[i])
                else if (val[m] - val[n] < -0.01 && step[m] - step[n] < -0.01 && haplotypes[index + m].offset[i] > haplotypes[index + n].offset[i])
                {
                    std::swap(haplotypes[index + m].offset[i], haplotypes[index + n].offset[i]);
                }
            }
        }
        
        for (int m = 0; m < ploidy; m++)
        {
            auto& newpriors = haplotypes[index + m].getnewprior(i);
            if (newpriors[0] == 0 || newpriors[1] == 0) continue;
            //step[m] = step[m] * (abssum / plainsum) + (step[m] - plainsum / ploidy) * (abssum / plainsum);             
            double centered = domaxcentered ? maxcentered : 2 * std::min(newpriors[0], newpriors[1]);
            // Was later:
            //if (burnin) centered = 1;
            
            //if (burnin) centered = 0;
            if (nocentertaper) centered = 1;
            if (simoffset && !burnin) centered = haplotypes[index + m].sim[i];
            if (mulsimoffset && !burnin) centered = haplotypes[index + m].sim[i];
            if (!mulsimoffset && crosssimoffset && !burnin) centered = 0;
            if (!tensionoffset) step[m] += centered * haplotypes[index + m].offset[i];
            double origstep = step[m];

            if (!snowball && !bigsnow)
            {
                haplotypes[index + m].momentum[i] *= momspeed;
                haplotypes[index + m].momentum[i] += step[m];
                step[m] = haplotypes[index + m].momentum[i];
            }
            else
            {
                step[m] += haplotypes[index + m].momentum[i];                
            }

            if (bigsnow)
            {
                haplotypes[index + m].momentum[i] = step[m] - origstep * (1 - momspeed);
                haplotypes[index + m].momentum[i] *= momspeed2;
            }
            
            step[m] = modifiedclamp ? std::clamp(step[m], -clamplim / stepsize, clamplim / stepsize) : std::clamp(step[m], -clamplim, clamplim);
            if (clampmom)
            {
                haplotypes[index + m].momentum[i] = modifiedclamp ? std::clamp<float>(haplotypes[index + m].momentum[i], -clamplim / stepsize, clamplim / stepsize) : std::clamp<float>(haplotypes[index + m].momentum[i], momstepclamp ? -newpriors[0] : -clamplim, momstepclamp ? newpriors[1] : clamplim);
            }

            if (snowball)
            {
                haplotypes[index + m].momentum[i] = step[m];
                haplotypes[index + m].momentum[i] *= momspeed;
            }
            
            if (noburninmom && burnin) haplotypes[index + m].momentum[i] = 0;

            val[m] += (step[m]) /* * (reads[i][0] + reads[i][1])*/ * stepsize * (noanypriorweight ? 1.0 : haplotypes[index + m].getanyprior(i)) * (1 + (stepszoffs ? haplotypes[index + m].offset[i] : 0)); // TODO: Needs to handle non-read count as well
            for (int j = 0; j < 2; j++)
            {
                float old = newpriors[j];
                int jstep = j ? -1 : 1;
                newpriors[j] = std::clamp<double>(exp(val[m] * jstep) / (exp(val[m] * jstep) + 1.0), updeps, 1 - updeps);
                double pseudostep = newpriors[j] - old;
                double otherstep = (fabs(j - ratio[m]) - old);
                if (minimalstep && pseudostep * step[m] * jstep <= 0 && (!guardminimum || otherstep * step[m] * jstep >= 0))
                {
                    //newpriors[j] = std::clamp<double>(nexttoward(old, old + step[m] * jstep), 1e-5, 1 - 1e-5);
                    //newpriors[j] = std::clamp<double>(nexttoward(old, old + jstep * (-1 * std::signbit(step[m]) + 1 * std::signbit(-step[m]))), 1e-5, 1 - 1e-5);
                    newpriors[j] = std::clamp<double>(nexttoward(old, old + jstep * (-1 * std::signbit(step[m]) + 1 * std::signbit(-step[m])) + (-1 * std::signbit(otherstep) + 1 * std::signbit(-otherstep))), updeps, 1 - updeps);
                }
            }
        }
    };

//#pragma omp parallel for schedule(dynamic, 120), num_threads(ploidy * 2)
#pragma omp taskloop num_tasks(ploidy * 2)
    for (int i = 0; i < genotypes.size(); i++)
    {
        if (updallpriors || reads[i][0] + reads[i][1] > 0 || genotypes[i] >= 0)
        {
            array<ratiotype, ploidy> ratio;
            double means[2] = {0};
            double binomsimcomppower = 1;
            if (liftmeannprior)
            for (int m = 0; m < ploidy; m++)
            {
                auto& priors = haplotypes[index + m].getprior(i);
                for (int j = 0; j < 2; j++)
                {
                    means[j] += priors[j];
                }                            
            }

            if (binomsimcomp)
            {
                double avsim = 0;
                for (int m = 0; m < ploidy; m++)
                {
                    for (int j = 0; j < ploidy; j++)
                    {
                        if (m == j) continue;
                        avsim += haplotypes[index + m].crosssim[i][j];
                    }
                }
                avsim /= (ploidy - 1) * ploidy;
                binomsimcomppower = 1 - avsim;
            }
            for (int m = 0; m < ploidy; m++)
            {
                auto reads = this->reads[i];
                double ourposteriormix = simposteriormix ? fabs(antisimposterior - haplotypes[index + m].sim[i]) :
                                        ((postmixred ? fabs(antisimposterior - haplotypes[index + m].sim[i]) : 1) * posteriormix);
                auto& priors = haplotypes[index + m].getprior(i);                                        

                if (newpostmix) ourposteriormix = 1 + haplotypes[index + m].sim[i] * (antiredcert ? certterm + priors[0] * priors[1] * certfactor : 1) * (-1 + posteriormix * (simredcert ? certterm + priors[0] * priors[1] * certfactor : 1));

                if (redcertmix)
                {
                    ourposteriormix *= certterm + priors[0] * priors[1] * certfactor;
                }

                array<ratiotype, ploidy> data[2];
                data[1].fill(0.f);

                bool now = true;
                data[1][0] = 1.0f;
                array<array<double, 2>, ploidy> mixposteriors;
                for (int j = 0; j < ploidy; j++)
                {
                    for (int n = 0; n < 2; n++)
                    {
                        double val;
                        if (arimeanmix) val = haplotypes[index + j].posteriorwo[i][n] * (1.0 - ourposteriormix) + haplotypes[index + j].posterior[i][n] * ourposteriormix;
                        else if (logitmeanmix)
                        {
                            val = pow(haplotypes[index + j].posteriorwo[i][n] / (haplotypes[index + j].posteriorwo[i][!n] + 1e-30f), 1.0 - ourposteriormix) *
                            pow(haplotypes[index + j].posterior[i][n] / (haplotypes[index + j].posterior[i][!n] + 1e-30f), ourposteriormix);
                            val = val / (1 + val);
                        }
                        else
                            val = pow(haplotypes[index + j].posteriorwo[i][n], 1.0 - ourposteriormix) * pow(haplotypes[index + j].posterior[i][n], ourposteriormix);

                        if (clampmix)
                        {
                            val = std::clamp<double>(val, updeps, 1 - updeps);
                        }

                        mixposteriors[j][n] = val;
                    }
                }

                for (int j = 0; j < ploidy; j++)
                {
                    if (j == m)
                    {
                        continue;
                    }             
                    now = !now;
                    data[now].fill(0.f);

                    auto& priors = haplotypes[index + j].getprior(i);
                    for (int k = 0; k < ploidy; k++)
                    {
                        //float sum = haplotypes[index + j].posterior[i][0] + haplotypes[index + j].posterior[i][1];
                        for (int n = 0; n < 2 && k + n < ploidy; n++)
                        {
                            data[now][k + n] += data[!now][k] * (burnin ? (mulpriorrest ? priors[n] : 0.5f) :
                            ((restposteriorwo ? 
                                mixposteriors[j][n] :
                            haplotypes[index + j].posterior[i][n]) * (antipriorpp ? priors[!n] : 1.0f))) /*/ /* sum*/;
                        }
                    }
                }

                double sums[2] = {0};
                for (int j = 0; j < 2; j++)
                {
                    for (int a = 0; a < ploidy; a++)
                    {
                        double base = data[now][a];
                        int counts[2] = {ploidy - 1 - a, a};
                        counts[j]++;
                        if (genotypes[i] != -1 && counts[1] != genotypes[i]) base *= std::max(domarkeps ? ourmap.otherepses[i] : 0.0f, epsothergeno) * (weakeneps ? (std::min(priors[j], priors[!j])) * 2 : 1.0f);

                        for (int k = 0; k < 2; k++)
                        {
                            if (reads[k] && !counts[k])
                            {
                                base *= 0;
                                continue;
                            }

                            if (!disableplacement)
                            {
                            int opts = reads[k] + counts[k] - 1;
                            // counts[k] groups, counts[k] -1 sentinel elements identifying borders
                            for (int z = 0; z < counts[k] - 1; z++)
                            {
                                base *= opts - z;
                                base /= counts[k] - 1 - z;
                            }
                        }
                        }
                        if ((burnin || permpostburnin) && !disableperm)
                        {
                            if (allhets)
                            {
                                if (counts[0] != ploidy || counts[1] != ploidy)
                                {
                                    base *= ploidy;
                                }
                            }
                            else
                            if (!altperm)
                            for (int z = 0; z < counts[0]; z++)
                            {
                                base *= ploidy - z;
                                base /= counts[0] - z;
                            }
                            else
                            {
                                int val = counts[j];
                                if (counts[j] < ploidy && ploidy - counts[j] < val) val = ploidy - counts[j];
                                base /= val;
                            }
                        }
                        int readsum = reads[0] + reads[1];
                        for (int k = 0; k < 2; k++)
                        {
                            for (int j = 0; j < reads[k]; j++)
                            {
                                base *= pow(counts[k], binomsimcomppower);
                                base /= ploidy * 0.5;
                                // Only one side of symmetry
                                if (!k)
                                {
                                    base *= readsum - j;
                                    base /= j + 1;
                    }
                            }
                        }
                        if ((burnin && mulpriorself) || (!burnin && propriorself)) base *= priors[j];
                        if (!burnin && antipriorself) base *= priors[!j];
                        if (!burnin && earlypost) base *= (selfposteriorwo && !nonwoearly) ?
                            (mixposteriors[m][j]) : haplotypes[index + m].posterior[i][j];
                        if ((burninassgn && burnin) || (postassgn && !burnin))
                        {
                            double notme = allnotme ? pow((ploidy - 1.0) / ploidy, reads[0] + reads[1]) :
                                                      (reads[j] ? pow((counts[j] - 1.0) / counts[j], reads[j]) : 1.0);
                            sums[j] += base * (1.0 - notme);
                            sums[0] += base * notme * ((earlypost && !burnin) ? ((selfposteriorwo && !selfpriorunass) ? haplotypes[index + m].posteriorwo[i][0] : haplotypes[index + m].posterior[i][0]) : 0.5);
                            sums[1] += base * notme * ((earlypost && !burnin) ? ((selfposteriorwo && !selfpriorunass) ? haplotypes[index + m].posteriorwo[i][1] : haplotypes[index + m].posterior[i][1]) : 0.5);
                        }
                        else
                        {
                            sums[j] += base;
                        }                        
                    } 
                    if (!burnin && !earlypost) sums[j] *= selfposteriorwo ? haplotypes[index + m].posteriorwo[i][j] : haplotypes[index + m].posterior[i][j];
                    if (liftmeannprior) sums[j] *= means[j] + (antiselfmean ? -priors[j] + 0.5 : 0.0);
                    //if (!burnin) sums[j] *= haplotypes[index + m].posterior[i][j] * priors[!j];
                    /*else
                    {
                        double factor = 0;
                        for (int k = 0; k < ploidy; k++)
                        {
                            auto& subpriors = haplotypes[index + k].getprior(i);
                            factor += subpriors[j];
                        }
                        factor /= ploidy;
                        sums[j] *= factor;
                    }*/
                    /*else
                        sums[j] *= priors[j];*/
                    /*else 
                        sums[j] * (reads[j] + 0.5) / (reads[0] + reads[1] + 1);*/
                }

                if (index == 0 && m == 0 && i == 11)
                {
                    printf("\n DATA: %lf %lf\t%f %f %f %f\t%f %f\n", sums[0], sums[1], data[now][0], data[now][1], data[now][2], data[now][3], haplotypes[index + 0].posterior[i][0], haplotypes[index + 0].posterior[i][1]);
                }
                ratio[m] = sums[0] / (sums[0] + sums[1] + 1e-30f);
            }
            updatenewpriors(i, ratio);
            /*if (!burnin)
            {
            for (int m = 0; m < ploidy; m++)
            {
                ratio[m] = haplotypes[index + m].posterior[i][0];
            }
            updatenewpriors(i, ratio);
            updatenewpriors(i, ratio);
            }*/
        }
    }
#pragma omp taskwait
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
    if (dovar || wounc) newanypriors = anypriors;
}

extern std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor> fwbw[ploidy];

/*#pragma omp threadprivate(fwbw)
std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor> fwbw[ploidy];*/

void doit()
{
    int hapnum = basehaps;    

    #pragma omp parallel
#pragma omp single
{
    std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor> fwbw[ploidy];
    #pragma omp taskloop num_tasks(24), private(hapnum, fwbw)
    //#pragma omp parallel for /*num_threads(16),*/ private(hapnum, fwbw)
    for (int i = 0; i < inds.size(); i++)
    {
        //#pragma omp parallel for num_threads(ploidy * 2), collapse(2), private(hapnum)
        //std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor>* fwbw = ::fwbw;

        hapnum = basehaps + i * ploidy;
        for (int k = 0; k < ploidy; k++)
        {
            haplotypes[hapnum + k].fwbw = &fwbw[k][0];
        }
        for (int k = 0; k < ploidy; k++)
        {
            for (int fw = 0; fw < 2; fw++)
            #pragma omp task firstprivate(i, k, fw, hapnum)
            {                            
                individ& ind = inds[i];
                haplotypes[hapnum + k].fwbw[fw].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                if (fw && fullwo)
                    haplotypes[hapnum + k].fwbw[2].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                if (fw /*&& ibdfactors*/)
                    haplotypes[hapnum + k].fwbw[2 + fullwo].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                if (fw /*&& ibdfactors*/ && nonsimfactor)
                    haplotypes[hapnum + k].fwbw[2 + fullwo + nonsimfactor].resize(haplotypes.size() * 2, ourmap.chromposes.size());    
                if (!burnin) haplotypes[hapnum + k].dofwbw(fw, ourmap);
            }
        }
        
        #pragma omp taskwait
        individ& ind = inds[i];
        bool flipped = !burnin && ind.handleflip(hapnum);
        if (!flipped)
        {
            //printf("Nudge %d/%d\n", hapnum, haplotypes.size());
            if (!burnin) ind.doposteriorhaplotypes(hapnum);
            ind.nudgehaplotypes(hapnum);
        }
        for (int k = 0; k < ploidy; k++)
        {
            haplotypes[hapnum + k].fwbw = nullptr;
        }
    }
    printf("GONNA WAIT FOR THE BITTER END\n");
    fflush(stdout);
    #pragma omp taskwait
}

    priors = newpriors;
    if (dovar || wounc) anypriors = newanypriors;
}

void readdummy(const char* mapname, const char* genoname, const char* allowedrefsname)
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
        prev += epsiloncM;
        if (pos < prev)
        {
            pos = prev;
        }
        ourmap.chromposes.push_back(pos);
        prev = pos;
    }
    ourmap.chromstarts.push_back(d);

    FILE* indfile = fopen(genoname, "rt");
    FILE* allowedrefsfile = fopen(allowedrefsname, "rt");
    int n;
    fscanf(indfile, "%d", &n);
    inds.resize(n);
    fscanf(allowedrefsfile, "%d", &n);
    for (individ& ind : inds)
    {
        ind.genotypes.resize(d);
        ind.reads.resize(d);
        std::fill(ind.reads.begin(), ind.reads.end(), std::array<int, 2>{0, 0});
        std::fill(ind.genotypes.begin(), ind.genotypes.end(), -1);
        for (int i = 0; i < d; i++)
        {            
            char tmp[255];
            fscanf(indfile, "%s", tmp);
            int a, b;
            if (sscanf(tmp, "%d;%d", &a, &b) == 2 && a >= 0 && b >= 0)
            {
                ind.reads[i]= {a, b};
            }
            else if (sscanf(tmp, "%d", &a) == 1)
            {
                ind.genotypes[i] = a;
            }
        }
        for (int i = 0; i < 8; i++)
        {
            int num;
            fscanf(allowedrefsfile, "%d", &num);
            ind.allowedrefs[i] = num;
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
                now.getprior(j)[val] = 1.f - refeps;
                now.getprior(j)[!val] = refeps;
                now.getanyprior(j) = true;
            }
            else
            {
                now.getanyprior(j) = false;
            }
        }
    }
}

void readerrors(const char* errorsname)
{
    FILE* errorsfile = fopen(errorsname, "rt");
    int n;
    fscanf(errorsfile, "%d", &n);
    ourmap.otherepses.resize(n);
    for (int i = 0; i < n; i++)
    {
        float val;
        fscanf(errorsfile, "%f", &val);
        if (val < 1)
        {
            val = std::min(1.0f, val / (ploidy - 1) * (1 / (1 - val)));
        }
        if (val < 0) val = 0;
        ourmap.otherepses[i] = val; 
    }
}

int main(int argc, char** argv) 
{/*
#ifdef _OPENMP
    omp_set_max_active_levels(2);
#endif*/
    //readdummy("potato_chr1.map", "potato_missing.gen");
    readdummy("polypop_1.map", "polyphref_1_0.gen", "allowedrefs__0.out");
    readerrors("polyerrors.out");
    readrefs("polyref_1_0.hap");
    int indi = 0;
    initinds();
    for (auto& ind : inds)
    {
        for (int i = 0; i < ind.genotypes.size(); i++)
        {
            if (ind.genotypes[i] == -1) continue;

            int mins[2] = {2, 2};
            int maxs[2] = {0, 0};
            for (int j = 0; j < basehaps; j++)
            {
                int index = j / (basehaps / 2);
                int val = haplotypes[j].getprior(i)[1] > 0.5;
                mins[index] -= 1 - val; //= std::min(val * 2, mins[index]);
                maxs[index] += val; //= std::max(val * 2, maxs[index]);
                mins[index] = std::max(0, mins[index]);
                maxs[index] = std::min(ploidy, maxs[index]);
            }
            mins[0] += mins[1];
            maxs[0] += maxs[1];

            if (ind.genotypes[i] < mins[0] || ind.genotypes[i] > maxs[0])
            {
                printf("NOT WORKING %d:%d %d (%d, %d)", indi, i,ind.genotypes[i], mins[0], maxs[0]);
                for (int j = 0; j < basehaps; j++)
                {
                    printf(" %d", haplotypes[j].getprior(i)[1] > 0.5);
                }
                printf("\n");
            }
        }
        indi++;
    }
    //inds.resize(2);
    double origstepsize = stepsize;
    burnin = true;
    stepsize = 0.2;
    for (int iter = 0; iter < 5000; iter++)
    {
        if (iter == 500)
    {
            burnin = false;
            stepsize = origstepsize;
        }
        for (int i = 0; i < inds.size(); i += inds.size() - 1)
        {
            for (int j = 0; j < 15; j++)
            {
                printf("%d %d", i, j);
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t%.3f %.3f", haplotypes[basehaps + i * ploidy + k].getprior(j)[1], haplotypes[basehaps + i * ploidy + k].getprior(j)[0]);
                }

                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t%.3f ", haplotypes[basehaps + i * ploidy + k].posterior[j][1]);
                }

                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t%.3f ", haplotypes[basehaps + i * ploidy + k].posteriorwo[j][1]);
                }

                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t% 01.3f ", haplotypes[basehaps + i * ploidy + k].offset[j]);
                    if (!burnin && iter >= 1750) haplotypes[basehaps + i * ploidy + k].offset[j] *= offsetdecay;
                }
                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t% 01.3f ", haplotypes[basehaps + i * ploidy + k].sim[j]);
                }

                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t%.3f ", 1 - haplotypes[basehaps + i * ploidy + k].desired[j]);
                }
                printf("\n");
            }
        }
        if (iter >= startstepshrink) stepsize *= 0.999;
        else
        if (!burnin && iter < endstepgrow) stepsize *= 1.0024;
        if (iter % tensionreset == 0) tension = 1.00;
        if (iter % tensionreset == 0 && postmixreset) posteriormix = 0;
        if (!burnin && !newNed)
        {
            Ne *= Nedecay;
            Ne -= Nestep;
            Ne = std::max<float>(Ne, Neend);
        }
        if (!burnin && newNed)
        {
            Ne = Neend + (Ne - Neend) * Nedecay;
        }
        tension *= tensiongrow;
        flipdrag += fldrstep;
        if (!burnin) posteriormix = 1.0 - (1.0 - posteriormix) * (1.0 - respostmix);
        printf("Test! %d %lf\n", iter, likelihood);
        likelihood = 0;
        doit();
        if (iter % 1000 == 999)
        {
    char filename[255];
    const char* letter = argv[1];

            sprintf(filename, "poly%s_%04d.vcflike", letter, iter);
    FILE* out = fopen(filename, "wt");
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

            sprintf(filename, "poly%s_%04d.out", letter, iter);
    out = fopen(filename, "wt");
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
    }
}
