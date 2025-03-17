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
const constexpr int maxreads = 20; 
const constexpr int permcount = ipow(ploidy, ploidy);
float stepsize = 0.05 / 3;
bool burnin = false;

template<class column> void doemit(column& c, float& anyprior, genprob& prior, int marker);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d);

vector<vector<genprob> > priors;
vector<vector<genprob> > newpriors;
vector<vector<float> > anypriors;

struct haplotype
{
    vector<genprob> posterior;
    vector<float> offset;

    ArrayXXf* fwbw;
    vector<double> renorm[2];
    genprob& getprior(int m) const;
    genprob& getnewprior(int m) const;
    float& getanyprior(int m) const;
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
                if (!fw /*&& getanyprior(from)*/) doemit(col, getanyprior(from), getprior(from), from);
                dotransition(col, col, themap, from, step);
            }
            if (fw /*&& getanyprior(m + sidestep)*/) doemit(col, getanyprior(m + sidestep), getprior(m + sidestep), m + sidestep);

            for (int i = 0, j = getindex() / ploidy * ploidy; i < ploidy; i++, j++)
            {
                col(j * 2) = 0;
                col(j * 2 + 1) = 0;
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

float& haplotype::getanyprior(int m) const
{
    return anypriors[m][getindex()];
}

int basehaps;

template<class column> void doemit(column& c, float& anyprior, genprob& prior, int marker)
{
    vector<genprob>& ourPrior = priors[marker];
    vector<float>& ourAnyPrior = anypriors[marker];
    #pragma ivdep
    for (int i = 0; i < haplotypes.size(); i++)
    {
        float old = c[i * 2] + c[i * 2 + 1];
        for (int j = 0; j < 2; j++)
        {
            float val = 0.f;
            val += (anyprior ? prior[j] : 1.0f) * ourPrior[i][j];
            
            float anyPriorW = /*anyprior * */ourAnyPrior[i] ? 1.0f : 0.0f;
            val *= anyPriorW;
            val += 0.5f * (1.0f - anyPriorW) * (anyprior ? prior[j] : 1.0f);
            c[i * 2 + j] = old * val;
//        if (val < 0 || val > 1) printf("%f\n", val);
        }
    }
}

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d)
{
    // Careful! c and c2 might coincide
    float dist = (themap.chromposes[marker + d] - themap.chromposes[marker]) * d * -0.02 * Ne;
    float nonrec = expf(dist);
    float rec = -expm1f(dist) / haplotypes.size();
    float sum = c.sum();
    float subsum;
    int prevbase = -1;
    for (int i = 0; i < haplotypes.size(); i++)
    {
        int base = i / ploidy * ploidy;
        if (base != prevbase)
        {
            subsum = c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1)).sum();
            prevbase = base;
        }
        float old = c[i * 2] + c[i * 2 + 1];
        for (int j = 0; j < 2; j++)
        {
            c2[i * 2 + j] = (old + subsum * 0.1) * nonrec + sum * rec;            
        }
    }
}

std::mt19937 rng;

struct individ
{
    vector<int> genotypes;
    vector<array<int, 2>> reads;
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
    std::uniform_real_distribution<float> distribution(-0.1, 0.1);

    for (int j = 0; j < ploidy; j++)
    {        
        haplotypes[index + j].posterior.resize(genotypes.size());
        haplotypes[index + j].offset.resize(genotypes.size());
        bool first = true;
        for (int i = 0; i < genotypes.size(); i++)
        {
            double genotype = genotypes[i];
            bool any = false;
            if (genotype < 0)
            {
                int readsum = reads[i][0] + reads[i][1];
                if (readsum)
                {
                    genotype = (reads[i][1] + (ploidy - 1) * 0.5) / (readsum + ploidy - 1) * ploidy;
                    haplotypes[index + j].getanyprior(i) = 1.0f - powf(powf(0.5, 1.0f / ploidy), readsum);
                }
                any = reads[i][0] && reads[i][1];
            }
            else
            {
                haplotypes[index + j].getanyprior(i) = true;
                any = genotype >= 1 && genotype <= ploidy - 1;
            }
            if (genotype >= 0)
            {                
                haplotypes[index + j].offset[i] = distribution(rng);
                float val = std::clamp<float>((genotype / 1.0f / ploidy) * (1.0f + haplotypes[index + j].offset[i]), 1e-5f, 1 - 1e-5f);
                haplotypes[index + j].getprior(i)[0] = 1.0f - val;
                haplotypes[index + j].getprior(i)[1] = val;                
                if (first && any && (j == 0 || j == ploidy - 1))
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

    #pragma omp parallel for schedule(guided, 100)
    for (int m = 0; m < haplotypes[index].fwbw[0].cols(); m++)
    {
        // TODO PRECALC ACCEL PLOIDY > 2
        bool first = true;
        double firstthisscore = 0;
        for (int p = permcount - 1; p >= 0; p--)
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
                    double firstscore = sum;
                    for (int j = 0; j < ploidy; j++)
                    {
                        firstscore += haplotypes[index + j].renorm[1][m];   
                        firstscore += haplotypes[index + j].renorm[0][m];
                    }
                    #pragma omp atomic
                    likelihood += firstscore;
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
        for (int p = permcount - 1; p >= 0; p--)
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
        array<float, ploidy> offset;
        for (int i = 0; i < haplotypes[index].posterior.size(); i++)
        {
            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                prior[j] = haplotypes[index + j].getprior(i);
                posterior[j] = haplotypes[index + j].posterior[i];
                offset[j] = haplotypes[index + j].offset[i];
            }

            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                int permval = i > bestmarker ? perm[j] : j;
                haplotypes[index + j].getnewprior(i) = prior[permval];
                haplotypes[index + j].posterior[i] = posterior[permval];
                haplotypes[index + j].offset[i] = offset[permval];
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
                //float weight = priors[m][k][0] * priors[m][k][0] + priors[m][k][1] * priors[m][k][1];
                for (int z = 0; z < 2; z++)
                {
                    haplotypes[index + j].posterior[m][z] += probs(k * 2 + z);
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
    auto updatenewpriors = [this, index] (int i, array<float, ploidy>& ratio)
    {
        array<double, ploidy> val, step;
        double abssum = 0;
        double plainsum = 0;
        for (int m = 0; m < ploidy; m++)
        {
            auto& priors = haplotypes[index + m].getnewprior(i);

            float midpoint = ratio[m];
            double num = std::clamp<double>(priors[0], 1e-10, 1.);
            double denom = std::clamp<double>(priors[1], 1e-10, 1.);
            val[m] = log(num/denom);
            step[m] = 1.0 / (exp(val[m]) + 1) + midpoint - 1.0 + haplotypes[index + m].offset[i];
            if (index == 0 && i == 3) printf("\n %d %d %d %lf %lf\n", index, m, i, val[m], step[m]);

            abssum += fabs(step[m]);
            plainsum += step[m];
        }

        plainsum = fabs(plainsum) + 1e-10f;

        for (int m = 0; m < ploidy; m++)
        {
            for (int n = 0; n < ploidy; n++)
            {
                if (m == n) continue;
                if (val[m] < 0 && step[m] - haplotypes[index + m].offset[i] < 0 && val[n] > 0 && step[n] - haplotypes[index + n].offset[i] > 0 && haplotypes[index + m].offset[i] > haplotypes[index + n].offset[i])
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
            step[m] = std::clamp(step[m], -1.0, 1.0);

            val[m] += step[m] * stepsize /** (reads[i][0] + reads[i][1])*/; // TODO: Needs to handle non-read count as well
            for (int j = 0; j < 2; j++)
            {
                newpriors[j] = std::clamp(exp(val[m] * (j ? -1 : 1)) / (exp(val[m] * (j ? -1 : 1)) + 1.0), 1e-3, 1 - 1e-3);
                //if (j) newpriors[j] = 1.0 - newpriors[j];
            }
        }
    };

#pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < genotypes.size(); i++)
    {
        if (genotypes[i] != -1)
        {
            int genotype = genotypes[i];

            float probs[2] = {0.f};
            for (int m = 0; m < ploidy; m++)
            {
                for (int l = 0; l < 2; l++)
                probs[l] += haplotypes[index + m].getprior(i)[l];
            }

            float sum = probs[0] + probs[1] + 1e-30f;
            probs[0] /= sum;
            probs[1] /= sum;

            if (i < 10 && index < 4) printf(" %.3f/%.3f", probs[1]*ploidy, (float) genotype);

            for (int m = 0; m < ploidy; m++)
            {
                array<float, ploidy + 1> probs[2] = {{0.f}, {0.f}};
                probs[0].fill(0.f);
                probs[1].fill(0.f);
                probs[1][0] = 1.f;
                int now = 1;
                for (int j = 0; j < ploidy; j++)
                {
                    if (j == m)
                    {
                        continue;
                    }
                    auto& priors = haplotypes[index + m].getprior(i);
                    now = !now;
                    std::fill(probs[now].begin(), probs[now].end(), 0.f);
                    for (int k = 0; k < ploidy; k++)
                    {
                        //float sum = haplotypes[index + j].posterior[i][0] + haplotypes[index + j].posterior[i][1];
                        for (int n = 0; n < 2 && k + n < ploidy; n++)
                        {
                            probs[now][k + n] += probs[!now][k] * priors[n] /** haplotypes[index + j].posterior[i][n] / sum*/;
                        }
                    }
                }

                auto& priors = haplotypes[index + m].getprior(i);                
                float val1 = (genotype ? probs[now][genotype - 1] : 0.f) * priors[1];
                float val2 = probs[now][genotype] * priors[0];
                //float diff = (val1 - val2) / (val1 + val2);

                //updatenewpriors(m, i, val1, val2);
            }
        }
        else
        if (reads[i][0] + reads[i][1]> 0)
        {
            array<float, ploidy> ratio;
            for (int m = 0; m < ploidy; m++)
            {
                auto reads = this->reads[i];
                array<float, ploidy> data[2];
                data[1].fill(0.f);

                bool now = true;
                data[1][0] = 1.0f;
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
                            data[now][k + n] += data[!now][k] * (burnin ? priors[n] : haplotypes[index + j].posterior[i][n]) /*/ /* sum*/;
                        }
                    }
                }

                double sums[2] = {0};
                auto& priors = haplotypes[index + m].getprior(i);
                for (int j = 0; j < 2; j++)
                {
                    for (int a = 0; a < ploidy; a++)
                    {
                        double base = data[now][a];
                        int counts[2] = {ploidy - 1 - a, a};
                        counts[j]++;
                        for (int k = 0; k < 2; k++)
                        {
                            if (reads[k] && !counts[k])
                            {
                                base *= 0;
                                continue;
                            }

                            int opts = reads[k] + counts[k] - 1;
                            // counts[k] groups, counts[k] -1 sentinel elements identifying borders
                            for (int z = 0; z < counts[k] - 1; z++)
                            {
                                base *= opts - z;
                                base /= counts[k] - 1 - z;
                            }
                        }
                        sums[j] += base;
                    }
                    //if (!burnin) sums[j] *= priors[j];
                    if (!burnin) sums[j] *= haplotypes[index + m].posterior[i][j];                    
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
                haplotypes[hapnum + k].fwbw[fw].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                if (!burnin) haplotypes[hapnum + k].dofwbw(fw, ourmap);
            }
        }
        hapnum = basehaps + i * ploidy;
        individ& ind = inds[i];
        bool flipped = !burnin && ind.handleflip(hapnum);
        if (!flipped)
        {
            //printf("Nudge %d/%d\n", hapnum, haplotypes.size());
            if (!burnin) ind.doposteriorhaplotypes(hapnum);
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
    readdummy("potato_chr1.map", "potato_reads.gen");
    //inds.resize(2);
    initinds();
    for (int k = 0; k < 1000; k++)
    {
        burnin = k < 100;
        for (int i = 0; i < 2; i++)
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
                    printf("\t%.3f %.3f ", haplotypes[basehaps + i * ploidy + k].posterior[j][1], haplotypes[basehaps + i * ploidy + k].posterior[j][0]);
                }

                printf("\t");
                for (int k = 0; k < ploidy; k++)
                {
                    printf("\t% 01.3f ", haplotypes[basehaps + i * ploidy + k].offset[j]);
                    if (!burnin) haplotypes[basehaps + i * ploidy + k].offset[j] *= 0.995;
                }
                printf("\n");
            }
        }
        if (!burnin) stepsize *= 1.002;
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
