#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include <numeric>
#include <tuple>
#include <math.h>
#include <omp.h>

#include "imputato_data.h"


// Borrowed https://stackoverflow.com/questions/17719674/c11-fast-constexpr-integer-powers
constexpr int64_t ipow_(int base, int exp){
  return exp > 1 ? ipow_(base, (exp>>1) + (exp&1)) * ipow_(base, exp>>1) : base;
}
constexpr int64_t ipow(int base, int exp){
  return exp < 1 ? 1 : ipow_(base, exp);
}

const constexpr int permcount = ipow(ploidy, ploidy);

using ArrayXPf = Eigen::Array<float, Eigen::Dynamic, ploidy>;

float stepsize = 0.0050;
bool burnin = false;

template<class column> void doemit(column& c, float& anyprior, genprob& prior, int marker, int* indices);

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d, int index, int majorclass);

int basehaps;

void haplotype::dofwbw(bool fw, const map& themap, int majorclass)
{
    ArrayXXf& myfwbw = fwbw[majorclass][fw];
    ArrayXf tempcol;
    tempcol.resize(myfwbw.rows());
    int colcount = myfwbw.cols();

    int start = fw ? 0 : colcount - 1;
    int end = fw ? myfwbw.cols() : 0;
    int step = fw ? 1 : -1;
    int sidestep = fw ? 0 : -1;

    if (fw || true)
    {
        auto col = myfwbw.col(start);
        for (int k = 0; k < myfwbw.rows(); k++)
        {
            col(k) = initclassweights[majorclass][classes[k / 2]];
        }
    }
    else
        myfwbw.col(start).fill(1.0f / myfwbw.rows());
    renorm[majorclass][fw][start] = 0.0f;

    if (filterrefs && allowedrefs && onlyref)
    {
        int half = ((getindex() - basehaps) % ploidy) >= (ploidy / 2);
        int halffactor = ploidy;

        if (!halfparinit || !fw)
        {
            half = 0;
            halffactor = ploidy * 2;
        }

        for (int i = 0; i < myfwbw.rows(); i++)
        {
            bool ok = false;
            for (int j = (half * (halffactor)); j < ((half + 1) * (halffactor)); j++)
            {
                if ((*allowedrefs)[j] == i / 2) ok = true;
            }

            if (!ok) myfwbw.col(start)(i) = 0;
        }
    }

    int indices[myfwbw.rows() / 2 / ibdgroupsize]; //ibd2sort

    for (int m = start; m != end; m += step)
    {
        auto col = myfwbw.col(m + sidestep);
        double srcrenorm = 0;            

        auto donz = [&] ()
        {
            int from = m - sidestep - 1;
            if (prelatenz && stepfrom) from -= sidestep;
            if (from >= 0 && from < myfwbw.cols() && (donzcmin || donzcess))
            {
                int nzc = 0;
                float nzmin = 1e30f;
                double nzsum = 0;
                double tempcolsum = 0;
                double sqsum = 0;
                for (int i = 0; i < myfwbw.rows(); i++)
                {
                    if (col(i))
                    {
                        nzmin = std::min(nzmin, col(i));
                        nzc++;
                    }
                    nzsum += col(i);
                    tempcolsum += tempcol(i);
                    if (donzcess)
                    {
                        sqsum += col(i) * col(i);
                    }
                }

                if (nzsum)
                {
                    if (donzcmin)
                    {
                        nzsum /= nzc;                        
                        nzmin *= nzminfactor * (1.0 - (nzsum - nzmin) / nzsum);
                        float scale = nzsum / (nzsum - nzmin);
                        for (int i = 0; i < myfwbw.rows(); i++)
                        {
                            if (col(i))
                            {
                                col(i) -= nzmin;
                                col(i) *= scale;
                            }
                        }
                    }
                    if (donzcess && nzc > 1)
                    {
                        double sim = 1;
                        if (sqsum)
                        {
                            double ess = nzsum * nzsum / sqsum;
                            sim = (ess - 1) / (nzc - 1);
                        }
                        sim *= sim;
                        sim *= nzminfactor;
                        if (sim < 1.0 - nzmaxfactor) sim = 1.0 - nzmaxfactor;
                        if (!isfinite(sim) || sim < 0 || sim > 1)
                        {
                            printf("NZWARN: %d %d %d %lf %lf %lf %d\n", getindex(), m, step, sim, nzsum, sqsum, nzc);
                        }
                        col = col * (1 - sim) + tempcol * sim * nzsum / tempcolsum;
                    }
                }
            }
        };
        if (m - sidestep - 1 >= 0)
        {
            int from = m - sidestep - 1;
            srcrenorm = renorm[majorclass][fw][from];
            myfwbw.col(m + sidestep) = myfwbw.col(from);
            if (!fw /*&& getanyprior(from)*/) doemit(col, getanyprior(from), getprior(from), from, indices);
            tempcol = 0.5 * (col + col.reshaped(2, col.size() / 2).colwise().reverse().reshaped());
            if (!fw && prelatenz) donz();            
            dotransition(col, col, themap, from, step, getindex(), majorclass);

            if (!latenz && (fw || !prelatenz)) donz();
        }

        for (int i = 0, j = getindex() / ploidy * ploidy; i < ploidy; i++, j++)
        {
            col(j * 2) = 0;
            col(j * 2 + 1) = 0;
        }

        if (fw /*&& getanyprior(m + sidestep)*/)
        {
            if (fullwo)
            {
                fwbw[majorclass][2].col(m + sidestep) = myfwbw.col(m + sidestep);
            }
            if (latenz) doemit(tempcol, getanyprior(m + sidestep), getprior(m + sidestep), m + sidestep, indices);            
            doemit(col, getanyprior(m + sidestep), getprior(m + sidestep), m + sidestep, indices);
        }

        if (latenz && (fw || !prelatenz)) donz();

        float sum = col.sum();
        sum += 1e-32;
        
        renorm[majorclass][fw][m + sidestep] = srcrenorm + log(sum);
        col *= expf(srcrenorm - renorm[majorclass][fw][m + sidestep]);
        if (fw && fullwo)
        {
            fwbw[majorclass][2].col(m + sidestep) *= expf(srcrenorm - renorm[majorclass][fw][m + sidestep]);
        }
    }
}

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

template<bool antigroup = ::antigroup, bool alreadyfactor = ibdfactors, size_t N = 1, class T, class T2>
void sortibd2b(ArrayXf& probs, std::array<ArrayXXf*, N> first, std::array<ArrayXXf*, N> second, int m, int* indices, ArrayXf& ysums, T&& factors, T2&& relevelfactors)
{
    double minlevelibd = ::minlevelibd * (mulminibd ? N : 1);
    double sum = 0;
    float max = 0;
    int groupcount = (haplotypes.size() - basehaps) / ibdgroupsize;
    if (ysums.size() < groupcount) ysums.resize(groupcount);
    if (factors.size() < groupcount) factors.resize(groupcount);
    probs.resize(first[0]->col(m).size());
    for (int x = 0; x < basehaps * 2; x++)
    {
        probs[x] = first[0]->col(m)[x] * second[0]->col(m)[x] * (relevel ? relevelfactors[x] : 1.0f);
        //if (relevelfactors[x] < 0.9 || relevelfactors[x] > 1.1 || !isfinite(relevelfactors[x])) printf("INCORRECT RELEVEL %f %d %d\n", relevelfactors[x], x, m);
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

template<class column> void dotransition(column& c, column& c2, const map& themap, int marker, int d, int index, int majorclass)
{
    // Careful! c and c2 might coincide
    float dist = (themap.chromposes[marker + d] - themap.chromposes[marker]) * d * -0.02 * Ne;
    dist = std::max(-maxexpdist, dist);
    float nonrec = expf(dist);
    float recbase = std::max(-expm1f(dist), 1e-5f);
    auto classweights = haplotypes[index].classweights[majorclass];
    if (d < 0)
    {
        for (int i = 0; i < numclasses; i++)
        {
            for (int j = 0; j < numclasses; j++)
            {
                classweights[i][j] = haplotypes[index].classweights[majorclass][j][i];
            }
        }
    }
    array<float, numclasses> sums;
    {
        array<double, numclasses> fullsums{0};
        for (int i = 0; i < haplotypes[index].classes.size(); i++)
        {
            double old = c[i * 2] + c[i * 2 + 1];
            fullsums[haplotypes[index].classes[i]] += old;
        }
        for (int i = 0; i < numclasses; i++)
        {
            sums[i] = fullsums[i];
        }
    }    
    float subsum = 0;
    float subunc = 0;
    float certf = 0;
    float partsum = 0;
    int prevbase = -1;    
    for (int i = 0; i < haplotypes.size(); i++)
    {
        bool selfref = onlyref && ((halfpar || halfparinit) && i < basehaps) && false;
        bool needsubsum = killibd2trans || selfref;
        if (needsubsum)
        {
            int base = i / ploidy * ploidy;
            if (base != prevbase)
            {
                subsum = c(Eigen::seq(base * 2, (base + ploidy) * 2 - 1)).sum();
                /*subsum = 0;
                for (int j = base * 2; j < (base + ploidy) * 2; j++)
                    subsum += c[j];*/
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
        const auto& myclassweights = classweights[haplotypes[index].classes[i]];

        float nowsum = 0;
        for (int j = 0; j < numclasses; j++)
        {
            nowsum += sums[j] * myclassweights[j];
        }
        for (int j = 0; j < 2; j++)
        {
            c2[i * 2 + j] = (1.0f - (filter ? filterlevel : 0)) * ((old * (certf + (1 - certf) * (killibd2trans ? (killibd2unc ? subunc : std::max(old, subsum - old) / (subsum + 1e-30f)) : 1 ))) * nonrec + nowsum * recbase);
            //c2[i * 2 + j] = (1.0f - (filter ? filterlevel : 0)) * (old  * nonrec + nowsum * rec);
        }
    }
}

std::mt19937 rng;


void individ::samplehaplotypes(int index)
{
    // Very crude, biased
    std::uniform_real_distribution<float> distribution(-offsmagn, offsmagn);

    for (int j = 0; j < ploidy; j++)
    {
        haplotypes[index + j].offset.resize(genotypes.size());
    }
    for (int i = 0; i < genotypes.size(); i++)
    {
        double sum = 0;
        int readsum = reads[i][0] + reads[i][1];
        float factor = (genotypes[i] != -1) ? 1 : 1 - pow((ploidy - 1.0) / ploidy, readsum);
        for (int j = 0; j < ploidy; j++)
        {
            haplotypes[index + j].offset[i] = distribution(rng) * factor;
            if (halfparinit && j % 2)
            {
                haplotypes[index + j].offset[i] = -haplotypes[index + j - 1].offset[i];
            }
            sum += haplotypes[index + j].offset[i];
        }
        sum /= ploidy;
        if (meanoffs)
        {
            for (int j = 0; j < ploidy; j++)
            {
                haplotypes[index + j].offset[i] -= sum;
                haplotypes[index + j].offset[i] *= ploidy / (ploidy - 1.0f);
            }
        }
    }
    

    for (int j = 0; j < ploidy; j++)
    {
        haplotypes[index + j].posterior.resize(genotypes.size());
        haplotypes[index + j].posteriorwo.resize(genotypes.size());        
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
                {
                    for (int k = 0; k < ploidy; k++)
                    {
                        
                    }
                    haplotypes[index + j].getanyprior(i) = false;
                }
                any = reads[i][0] && reads[i][1];
            }
            else
            {
                any = genotype >= 1 && genotype <= ploidy - 1;
                genotype = (genotype + 0.005) / (ploidy + 0.01) * ploidy;
                haplotypes[index + j].getanyprior(i) = true;                
            }

            haplotypes[index + j].momentum[i] = 0;

            if (genotype >= 0)
            {
                if (firstatall == -1) firstatall = i;
            
                bool halfed = sampoffshalf && genotype < ploidy / 2;
                if (halfed)
                {
                    genotype = ploidy - genotype;
                    haplotypes[index + j].offset[i] = -haplotypes[index + j].offset[i];
                }
                double val = std::clamp<double>((genotype / 1.0f / ploidy) * (1.0f - haplotypes[index + j].offset[i]), updeps, 1.0 - updeps);
                if (halfed)
                {
                    val = 1.0f - val;
                    genotype = ploidy - genotype;
                    haplotypes[index + j].offset[i] = -haplotypes[index + j].offset[i];
                }
                //haplotypes[index + j].offset[i] = distribution(rng);//-haplotypes[index + j].offset[i];
                haplotypes[index + j].getprior(i)[0] = 1.0 - val;
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
    vector<array<double, permcount> > onescores;
    scores.resize(haplotypes[index].fwbw[0][0].cols());
    if (oneflip) onescores.resize(haplotypes[index].fwbw[0][0].cols());
    int indices[haplotypes.size() / ibdgroupsize];
    ArrayXf probs[ploidy], oneprobs[ploidy], ysums, factors;
    individ& ind = inds[(index - basehaps) / ploidy];
//    #pragma omp parallel for schedule(guided, 100), private(indices, probs, ysums, factors), shared(scores, ind), num_threads(ploidy * 2)
    #pragma omp taskloop num_tasks(ploidy * 2), private(indices, probs, oneprobs, ysums, factors), shared(ind, scores, onescores)
    for (int m = 0; m < haplotypes[index].fwbw[0][0].cols(); m++)
    {
        bool first = true;
        double firstthisscore = 0;
        double firstagnscore = 0;
        double firstthisscoreone = 0;
        int firstpow2 = 0;
        int firstpow2one = 0;
        float sims[ploidy][ploidy];
        float corrs[ploidy];

        for (int majorclass = 0; majorclass < haplotypes[index].herenummajor; majorclass++)
        {
            if (ibdfactors)
            {
                if (relevel)
                    for (int k = 0; k < (multirelev ? ploidy : 1); k++)
                    {
                        haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m).fill(1.0f);
                    }
                for (int k = 0; k < (antifactors ? ploidy : 1); k++)
                {
                    constexpr int count = antifactors ? 1 : ploidy;
                    std::array<ArrayXXf*, count> first;
                    std::array<ArrayXXf*, count> firstwo;
                    std::array<ArrayXXf*, count> second;

                    for (int j = 0; j < count; j++)        
                    {                    
                        first[j] = &haplotypes[index + k + j].fwbw[majorclass][1 + powofactors];
                        firstwo[j] = &haplotypes[index + k + j].fwbw[majorclass][1 + fullwo];
                        second[j] = &haplotypes[index + k + j].fwbw[majorclass][0];
                    }
                    sortibd2b<antigroup, false, count>(probs[k], first, second, m, indices, ysums, haplotypes[index + k].fwbw[majorclass][2 + fullwo].col(m), haplotypes[index + (multirelev ? k : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                    if (nonsimfactor) sortibd2b<antigroup, false, count>(probs[k], firstwo, second, m, indices, ysums, haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor].col(m), haplotypes[index + (multirelev ? k : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                }

                if (relevel)
                {
                    float sums[ploidy];
                    for (int k = 0; k < ploidy; k++)
                    {
                        double sum = 0;
                        for (int i = 0; i < haplotypes.size() * 2; i++)
                        {
                            sum += probs[k][i];
                        }
                        sums[k] = 1.0 / std::max(sum, 1e-30);
                    }
                    for (int i = 0; i < haplotypes.size(); i++)
                    {
                        float val = 0;
                        float vals[ploidy];
                        for (int k = 0; k < ploidy; k++)
                        {
                            vals[k] = (probs[k][i * 2] + probs[k][i * 2 + 1]) * sums[k];
                            vals[k] = std::min(vals[k], 1.0f - 1e-4f);
                            val += vals[k];
                        }
                        float origval = val;
                        if (!ind.singlerelevel)
                        {
                            if (val < 1.0f) val = 1.0f;
                            else
                            {
                                if (val > 1.999f) val = 1.999f;
                                val = (val - 1) / (2 - val);
                                val = 1 / (1 + val);
                                
                            }
                        }
                        else
                        {
                            if (val < ind.singlerelevel)
                            {
                                val = 1.0f;
                                /*if (multirelev) // reset to 1 per default
                                {
                                    for (int k = 0; k < ploidy; k++)
                                    {
                                        haplotypes[index + k].fwbw[2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2) = val;
                                        haplotypes[index + k].fwbw[2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2 + 1) = val;
                                    }
                                }*/
                            }
                            else
                            {                            
                                val = 1.0f;
                                for (int k = 0; k < ploidy; k++)
                                {
                                    float clampval = origval;
                                    if (clampval > ind.singlerelevel + vals[k] * ind.singlerelevel)
                                    {
                                        clampval = ind.singlerelevel + vals[k] * ind.singlerelevel;
                                    }
                                    //float newval = (vals[k] - 1) * (vals[k] + singlerelevel - 1) / (vals[k] * (2 * vals[k] + singlerelevel - 1));
                                    // ax / (ax + 1 - x) = c - (b - x)
                                    float newval = (vals[k] - 1) * (clampval - ind.singlerelevel - vals[k]) / (vals[k] * (clampval - ind.singlerelevel - vals[k] + 1));
                                    if (newval < 1e-3f) newval = 1e-3f;
                                    if (!isfinite(newval)) newval = 1.0f;
                                    if (multirelev)
                                    {
                                        haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2) = newval;
                                        haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2 + 1) = newval;
                                        if ((newval < 1.0f && index == 68 && m == 68) || newval < 0) printf("RELEVEL %d %d %d %d %f %f\t%f %f %f %f\n", index, m, i, k, newval, origval, vals[0], vals[1], vals[2], vals[3]);
                                    }
                                    if (newval < val) val = newval;
                                }
                            }
                        }
                        if (!multirelev)
                        {
                            if ((val < 1.0f && index == 68 && m == 68) || val < 0) printf("RELEVEL %d %d %d %f %f\t%f %f %f %f\n", index, m, i, val, origval, vals[0], vals[1], vals[2], vals[3]);
                            haplotypes[index].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2) = val;
                            haplotypes[index].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2 + 1) = val;
                        }

                        // TODO: Integrate both types better
                        if (ind.otherrelevel && multirelev /*&& (ind.genotypes[m] < 1 || ind.genotypes[m] > ploidy - 1)*/ && ind.genotypes[m] == -1)
                        {
                            for (int k = 0; k < ploidy; k++)
                            {
                                if (vals[k])
                                {
                                    float newval = 1.0f + (origval - vals[k]) / (vals[k]) * ind.otherrelevel;
                                    if (newval < 1e-3f) newval = 1e-3f;
                                    //if (newval > ploidy) newval = ploidy;
                                    if (!isfinite(newval)) newval = 1.f;
                                    
                                    haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2) *= newval;
                                    haplotypes[index + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m)(i * 2 + 1) *= newval;
                                }
                            }
                        }
                    }
                }
            }
        }

        // m != 0 just to get stats that are only computed at m == 0
        if (!updallpriors && genotypes[m] == -1 && reads[m][0] + reads[m][1] == 0 && m != 0)
        {
            for (int p = permcount - 1; p >= 0; p--)
            {
                scores[m][p] = -1.1e30f;
                if (oneflip) onescores[m][p] = -1.1e30f;
            }
            continue;
        }

        //#pragma ivdep
        for (int j = 0; j < ploidy; j++)        
        {
            int majorclass = haplotypes[index + j].mainmajorclass;
            if (!expklsim)
            {
                if (sortsim && groupibd2)
                {
                    sortibd2b(probs[j], {&haplotypes[index + j].fwbw[majorclass][1 + simpowo]}, {&haplotypes[index + j].fwbw[majorclass][0]}, m, indices, ysums, haplotypes[index + (antifactors ? j : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor].col(m), haplotypes[index + (multirelev ? j : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                }
                else
                {
                    probs[j] = haplotypes[index + j].fwbw[majorclass][1 + simpowo].col(m) * haplotypes[index + j].fwbw[majorclass][0].col(m);
                    if (sortsim)
                    {
                        sortibd2(probs[j], indices);
                    }
                }
                //corrs[j] = 1.0f / (sqrt((probs[j] * probs[j]).sum()) + 1e-30f);
                corrs[j] = sqrt(1.0f / probs[j].sum());
                //corrs[j] = 1.0f / (probs[j].sum() + 1e-30f);
            }
            else
            {
                probs[j] = haplotypes[index + j].fwbw[majorclass][1 + simpowo].col(m) * haplotypes[index + j].fwbw[majorclass][0].col(m) + 1e-30f;
                if (sortsim) sortibd2(probs[j], indices);
                corrs[j] = 1.0f / probs[j].sum();
            }
        }

        ind.maxshared[m] = 0;
        float sums[ploidy];
        for (int k = 0; k < ploidy; k++)
        {
            double sum = 0;
            for (int i = 0; i < haplotypes.size() * 2; i++)
            {
                sum += probs[k][i];
            }
            sums[k] = 1.0 / std::max(sum, 1e-30);
        }
        for (int i = 0; i < haplotypes.size(); i++)
        {
            float fraction = 0;
            for (int j = 0; j < ploidy; j++)
            {
                fraction += (probs[j][i * 2] + probs[j][i * 2 + 1]) * sums[j];
            }            
            if (fraction > ind.maxshared[m])
            {
                ind.maxshared[m] = fraction;
                ind.maxsharedid[m] = i;
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
                    haplotypes[index + j].crosssim[m][k] = 1 - std::clamp<float>((probs[j] * probs[k]).sum() * corrs[j] * corrs[k] * corrs[j] * corrs[k], 0, 1);
                    continue;
                }
                if (!expklsim)
                {
                    haplotypes[index + j].crosssim[m][k] = std::clamp<float>(sqrt(probs[j] * probs[k]).sum() * corrs[j] * corrs[k], 0, 1);
                    //haplotypes[index + j].crosssim[m][k] = 1 - std::clamp<float>(abs(probs[j] * corrs[j] - probs[k] * corrs[k]).sum() * 0.5, 0, 1);
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
        int onepow2s[ploidy][ploidy];
        double onesinglescores[ploidy][ploidy];

        for (int j = 0; j < ploidy; j++)
        {
            int majorclass = haplotypes[index + j].mainmajorclass;
            for (int k = 0; k < ploidy; k++)
            {
                double sumterm = 0;
                double sumoneterm = 0;
                if (oneflip)
                {
                    auto col = haplotypes[index + j].fwbw[majorclass][2 + fullwo + nonsimfactor + 1].col(m);
                    haplotypes[index + j].fwbw[majorclass][2 + fullwo + nonsimfactor + 1].col(m) = haplotypes[index + j].fwbw[majorclass][2].col(m);
                    doemit(col, haplotypes[index + k].getanyprior(m), haplotypes[index + k].getprior(m), m, indices);
                }
                if (sortflip && groupibd2)
                {
                    if (!vetoflip && !crossibd)
                    {
                        if (oneflip)
                        {
                            sortibd2b<antigflip>(oneprobs[k], {&haplotypes[index + j].fwbw[majorclass][2 + fullwo + nonsimfactor + 1]}, {&haplotypes[index + j].fwbw[majorclass][0]}, m, indices, ysums, haplotypes[index + (antifactors ? j : 0)].fwbw[majorclass][2 + fullwo].col(m), haplotypes[index].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                        }

                        sortibd2b<antigflip>(probs[k], {&haplotypes[index + j].fwbw[majorclass][1]}, {&haplotypes[index + k].fwbw[majorclass][0]}, m, indices, ysums, haplotypes[index + (antifactors ? j : 0)].fwbw[majorclass][2 + fullwo].col(m), haplotypes[index].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                    }
                    else
                    {
                        if (oneflip)
                        {
                            sortibd2b<antigflip, false>(oneprobs[k], {&haplotypes[index + j].fwbw[majorclass][2 + fullwo + nonsimfactor + 1]}, {&haplotypes[index + j].fwbw[majorclass][0]}, m, indices, ysums, factors, haplotypes[index + (multirelev ? j : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                        }

                        sortibd2b<antigflip, false>(probs[k], {&haplotypes[index + j].fwbw[majorclass][1]}, {&haplotypes[index + k].fwbw[majorclass][0]}, m, indices, ysums, factors, haplotypes[index + (multirelev ? j : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                    }
                }
                else
                {
                    abort();
                    probs[k] = haplotypes[index + j].fwbw[majorclass][1].col(m) * haplotypes[index + k].fwbw[majorclass][0].col(m);
                    if (sortflip) sortibd2(probs[k], indices);
                }
                // TODO returnera summan så vi har den
                for (int i = 0; i < haplotypes.size() * 2; i++)
                {
                    /*double a = haplotypes[index + j].fwbw[1].col(m)(i);
                    double b = haplotypes[index + k].fwbw[0].col(m)(i);*/

                    sumterm += probs[k](i);
                    if (oneflip) sumoneterm += oneprobs[k](i);
                }
                if (vetoflip && j != k)
                {
                    for (int z : {j, k})
                    {
                        double sumterm2  = 0;
                        double sumoneterm2  = 0;
                        if (oneflip)
                        {
                            sortibd2b<antigflip>(oneprobs[k], {&haplotypes[index + j].fwbw[majorclass][2 + fullwo + nonsimfactor + 1]}, {&haplotypes[index + j].fwbw[majorclass][0]}, m, indices, ysums, haplotypes[index + (antifactors ? z : 0)].fwbw[majorclass][2 + fullwo].col(m), haplotypes[index + (multirelev ? z : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                        }
                        sortibd2b<antigflip>(probs[k], {&haplotypes[index + j].fwbw[majorclass][1]}, {&haplotypes[index + k].fwbw[majorclass][0]}, m, indices, ysums, haplotypes[index + (antifactors ? z : 0)].fwbw[majorclass][2 + fullwo].col(m), haplotypes[index + (multirelev ? z : 0)].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));

                        for (int i = 0; i < haplotypes.size() * 2; i++)
                        {
                            /*double a = haplotypes[index + j].fwbw[1].col(m)(i);
                            double b = haplotypes[index + k].fwbw[0].col(m)(i);*/

                            sumterm2 += probs[k](i);
                            if (oneflip) sumoneterm2 += oneprobs[k](i);
                        }
                        if (sumterm2 < sumterm) sumterm = sumterm2;
                        if (sumoneterm2 < sumoneterm) sumoneterm = sumoneterm2;
                    }
                }
                singlescores[j][k] = log(frexp(sumterm + 1e-300, &pow2s[j][k])) + (haplotypes[index + k].renorm[majorclass][0][m] - haplotypes[index + j].renorm[majorclass][0][m]);
                if (oneflip) onesinglescores[j][k] = log(frexp(sumoneterm + 1e-300, &onepow2s[j][k]));
            }
        }

        for (int p = permcount - 1; p >= 0; p--)
        {
            array<int, ploidy> perm;
            bool badperm = !getploidyperm(p, perm);
            if (!badperm && (halfpar || halfparinit) && false)
            {
                for (int i = 0; i < ploidy; i++)
                {
                    if ((perm[i] < ploidy / 2) ^ (i < ploidy / 2)) badperm = true;
                }
            }
            if (badperm)
            {
                scores[m][p] = -1.1e30f;
                if (oneflip) onescores[m][p] = -1.1e30f;
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
            double sumone = 0;
            int pow2 = 0;
            int pow2one = 0;
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
                                    terms[z] += haplotypes[nowindex].fwbw[haplotypes[nowindex].mainmajorclass][k].col(m)(subindex);
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
                    if (oneflip)
                    {
                        sumone += onesinglescores[j][perm[j]];
                        pow2one += onepow2s[j][perm[j]];                    
                    }
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
                        int majorclass = haplotypes[index + j].mainmajorclass;
                        firstscore += haplotypes[index + j].renorm[majorclass][1][m];   
                        firstscore += haplotypes[index + j].renorm[majorclass][0][m];
                    }
                    firstscore += log(2) * pow2;
                    #pragma omp atomic
                    likelihood += firstscore;
                    haplotypes[index].likelihood = firstscore;
                }
                sum += flipbias;
                sumone += flipbias;
                firstpow2 = pow2;
                firstpow2one = pow2one;
                firstthisscore = sum;
                firstthisscoreone = sumone;
                firstagnscore = sumagn;
                //if (fabs(firstthisscore - firstthisscoreone + log(2) * (firstpow2 - firstpow2one)) > 0.0001) printf("BEJ %4d %4d %4d %d %d %lf %lf\n", index, m, p, firstpow2, firstpow2one, firstthisscore, firstthisscoreone);
                first = false;
            }

            //if (index == 16) printf("Flip: %d %d %d %f\n", index, m, p, sum);
            sum += log(2) * (pow2 - firstpow2);
            if (oneflip) sumone += log(2) * (pow2one - firstpow2one);
            scores[m][p] = sum - firstthisscore;
            //if (perm[0] >= 2 || perm[1] >= 2) scores[m][p] = -1.1e30f;
            if (oneflip) onescores[m][p] = sumone - firstthisscoreone;
            if (bothagn && sumagn < firstagnscore) scores[m][p] = sumagn - firstagnscore;
        }
    }

    #pragma omp taskwait

    int bestmarker = 0;
    int bestp = 0;
    double bestscore = -1.1e30f;
    double realbestscore = -1.1e30f;    
    for (int m = 0; m < haplotypes[index].fwbw[haplotypes[index].mainmajorclass][0].cols(); m++)
    {
        for (int p = permcount - 1; p >= 0; p--)
        {
            for (int mode = 0; mode < 1 + oneflip; mode++)
            {
                double sum = mode ? onescores[m][p] : scores[m][p];
                sum = sum * (noiseflippos ? std::uniform_real_distribution<double>(0, 1)(rng) : 1);
                if (sum < -1e30f) continue;

                if (!randflippos || bestscore < 0)
                {
                    if (sum > bestscore)
                    {
                        bestscore = sum;
                        realbestscore = sum;
                        bestp = p + mode * 1048576;
                        bestmarker = m;
                    }
                }
                else if (sum >= 0)
                {
                    if (!std::bernoulli_distribution(1 / (exp(sum - bestscore) + 1))(rng))
                    {                
                        bestp = p + mode * 1048576;
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
    }    

    return {bestmarker, bestp, realbestscore};
}

bool individ::handleflip(int index)
{
    auto [bestmarker, bestp, bestscore] = findflip(index);

    array<int, ploidy> perm;
    bool straight = true;
    getploidyperm(bestp & 1048575, perm);
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
        if (bestp > 1048575)
        {
            printf(" (ONEFLIP)");
        }
        else
        {
            printf(" (NoNEFLIP)");
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
                    //anyprior[j] = haplotypes[index + j].getanyprior(i);
                    newanyprior[j] = haplotypes[index + j].getnewanyprior(i);
                }
            }

            #pragma ivdep
            for (int j = 0; j < ploidy; j++)
            {
                int permval = j;
                if (bestp > 1048575)
                {
                    if (i == bestmarker) permval = perm[j];
                }
                else
                {
                    if (i > bestmarker) permval = perm[j];
                }
                
                haplotypes[index + j].getnewprior(i) = prior[permval];
                haplotypes[index + j].posterior[i] = posterior[permval];
                haplotypes[index + j].posteriorwo[i] = posteriorwo[permval];
                haplotypes[index + j].offset[i] = offset[permval];
                haplotypes[index + j].sim[i] = sim[permval];
                haplotypes[index + j].desired[i] = desired[permval];
                haplotypes[index + j].momentum[i] = momentum[permval];                
                if (dovar || wounc)
                {
                    //haplotypes[index + j].getanyprior(i) = anyprior[permval];
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
        for (int m = 0; m < haplotypes[index].fwbw[haplotypes[index].mainmajorclass][0].cols(); m++)
        {
            int mainmajorclass = haplotypes[index + j].mainmajorclass;
            if (ibd2sort && groupibd2)
            {
                sortibd2b(probs, {&haplotypes[index + j].fwbw[mainmajorclass][1]}, {&haplotypes[index + j].fwbw[mainmajorclass][0]}, m, indices, ysums, haplotypes[index + (ibdfactors && !antifactors ? 0 : j)].fwbw[mainmajorclass][2 + fullwo].col(m), haplotypes[index + (multirelev ? j : 0)].fwbw[mainmajorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
            }
            else
            {
                probs = haplotypes[index + j].fwbw[mainmajorclass][1].col(m) * haplotypes[index + j].fwbw[mainmajorclass][0].col(m);
                if (ibd2sort)
                {
                    sortibd2(probs, indices);
                }
            }

            if (fullwo)
            {
                if (ibd2sort && groupibd2)
                {
                    sortibd2b(probswo, {&haplotypes[index + j].fwbw[mainmajorclass][2]}, {&haplotypes[index + j].fwbw[mainmajorclass][0]}, m, indices, ysums, haplotypes[index + (ibdfactors && !antifactors ? 0 : j)].fwbw[mainmajorclass][2 + fullwo + nonsimfactor].col(m), haplotypes[index + (multirelev ? j : 0)].fwbw[mainmajorclass][2 + fullwo + nonsimfactor + oneflip + relevel].col(m));
                }
                else
                {
                    probswo = haplotypes[index + j].fwbw[mainmajorclass][2].col(m) * haplotypes[index + j].fwbw[mainmajorclass][0].col(m);
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
                    haplotypes[index + j].getnewanyprior(m) = std::min(1.0, std::max(1e-30, dovarunc ? uncvar : (onlyvar ? var : (1 - var + uncvar))));
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
    #pragma omp taskwait
}

void individ::nudgehaplotypes(int index)
{
    auto updatenewpriors = [this, index] (int i, array<ratiotype, ploidy>& ratio, float speed = 1.0)
    {
        ratiotype renormproportion = -1;
        array<double, ploidy> val, step, midpoints, origmidpoints, powomidpoints;
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

        if (renormproportion >= 0)
        {
            ratiotype sum = 0;
            for (ratiotype r : ratio)
            {
                sum += r;
            }

            sum /= ploidy;
            if ((sum < updeps && renormproportion < updeps) || (sum > 1 - updeps && renormproportion > 1 - updeps))
            {            
            }
            else
            {
                renormproportion = std::clamp<ratiotype>(renormproportion, updeps, 1 - updeps);
                sum = std::clamp<ratiotype>(sum, updeps, 1 - updeps);
                renormproportion = 1 / (1 - renormproportion);
                sum = 1 / (1 - sum);
                renormproportion /= sum;
                for (ratiotype& r : ratio)
                {
                    r *= renormproportion;
                }
            }
        }

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
                double midpointsval =  log(std::clamp<double>(midpoint, 1e-10, 1.)) - log(std::clamp<double>(midpoint1m, 1e-10, 1.));
                if (priorpowo) origmidpoints[m] = midpointsval;
                if (bothpowo) powomidpoints[m] = midpointsval;
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
                auto& priorsm = haplotypes[index + m].getprior(i);
                for (int k = 0; k < ploidy; k++)
                {
                    if (m == k) continue;//break; //continue;
                    auto& newpriorsk = haplotypes[index + k].getnewprior(i);
                    auto& priorsk = haplotypes[index + k].getprior(i);
                    if (caponpriors)
                    {
                        if (newpriorsm[0] == 0 || newpriorsk[0] == 0 || newpriorsm[1] == 0 || newpriorsk[1] == 0) continue;
                        if (priorsm[0] == 0 || priorsk[0] == 0 || priorsm[1] == 0 || priorsk[1] == 0) continue;
                    }
                    else
                    {
                        if (fabs(origmidpoints[m]) > 30 || fabs(origmidpoints[k]) > 30) continue;
                    }

                    // 1 in numerator implies switched order
                    double diff = plaincsdiff ? 1/(1+exp(origmidpoints[k])) - 1/(1+exp(origmidpoints[m])) : (origmidpoints[m] - origmidpoints[k]);
                    if (selflipcs && (midpoints[m] - midpoints[k]) * diff < 0) diff = -diff;
                    if (selzerocs && (midpoints[m] - midpoints[k]) * diff < 0) diff = 0;

                    if (bothpowo)
                    {
                        double origdiff = diff;
                        diff = plaincsdiff ? 1/(1+exp(powomidpoints[k])) - 1/(1+exp(powomidpoints[m])) : (powomidpoints[m] - powomidpoints[k]);
                        if (selflipcs && (midpoints[m] - midpoints[k]) * diff < 0) diff = -diff;
                        if (selzerocs && (midpoints[m] - midpoints[k]) * diff < 0) diff = 0;

                        if (fabs(origdiff) > fabs(diff)) diff = origdiff;
                    }
                    if (rewcsdiff) diff = diff / (1.0001 - haplotypes[index + m].crosssim[i][k]);
                    double delta = (diff + (midpoints[m] - midpoints[k] < 0 ? -1 : 1) * csbump * (invcsbump ? (1.0001 - haplotypes[index + m].crosssim[i][k]) - 1.0 + 0.0001 : 1)) * (fabs(haplotypes[index + m].offset[i]) + fabs(haplotypes[index + k].offset[i]))
                        * (invcs ? 1.0 / (1.0001 - haplotypes[index + m].crosssim[i][k]) - 1.0 + 0.0001 : haplotypes[index + m].crosssim[i][k])
                        * csscale;
                    midpoints[m] += delta;
                    //midpoints[k] -= delta;                    
                }
               midpoints[m] = std::clamp<double>(midpoints[m], -midpointcap, midpointcap);
               //ratio[m] = exp(midpoints[m]) / (1 + exp(midpoints[m]));
               ratio[m] = 1 / (1 + exp(-midpoints[m]));
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
            float nowspeed = speed;
            if (!burnin) nowspeed *= (1 - haplotypes[index + m].crosssim[i][m]);
            if (!simplestep)
            {
                if (logitstep)
                {
                    step[m] = (midpoint - val[m]) * nowspeed;
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
                step[m] = (midpoint - val[m]) * nowspeed;
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
            if (mulsimoffset && !burnin) centered *= haplotypes[index + m].sim[i];
            if (!mulsimoffset && crosssimoffset && !burnin) centered = 0;
            if (!tensionoffset && !lateoffset) step[m] += centered * haplotypes[index + m].offset[i];
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

            if (!tensionoffset && lateoffset) step[m] += centered * haplotypes[index + m].offset[i];

            val[m] += (step[m]) /* * (reads[i][0] + reads[i][1])*/ * stepsize * (noanypriorweight ? 1.0 : haplotypes[index + m].getanyprior(i)) * (1 + (stepszoffs ? haplotypes[index + m].offset[i] : 0)); // TODO: Needs to handle non-read count as well
            for (int j = 0; j < 2; j++)
            {
                float old = newpriors[j];
                int jstep = j ? -1 : 1;
                newpriors[j] = std::clamp<double>(1.0 / (exp(-val[m] * jstep) + 1.0), updeps, 1 - updeps);
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
        if (!(updallpriors || reads[i][0] + reads[i][1] > 0 || genotypes[i] >= 0)) continue;

        const auto reads = this->reads[i];
        array<ratiotype, ploidy> ratio;
        double means[2] = {0};
        double binomsimcomppower = 1;

        array<ratiotype, ploidy + 1> data[2];
        array<ratiotype, ploidy + 1> datatarget[2];
        bool now = true;        
        if (!burnin)
        {
            data[now].fill(0.f);
            data[now][0] = 1.0f;
            datatarget[now].fill(0.f);
            datatarget[now][0] = 1.0f;
            for (int j = 0; j < ploidy; j++)
            {
                now = !now;
                data[now].fill(0.f);
                datatarget[now].fill(0.f);

                auto& priors = haplotypes[index + j].getprior(i);
                for (int k = 0; k < ploidy; k++)
                {
                    //float sum = haplotypes[index + j].posterior[i][0] + haplotypes[index + j].posterior[i][1];
                    for (int n = 0; n < 2 && k + n <= ploidy; n++)
                    {
                        float sim = haplotypes[index + j].sim[i];
                        data[now][k + n] += data[!now][k] * /*pow(haplotypes[index + j].posteriorwo[i][n], 1 - sim) /** globallelebias[i][n]*/ /* *
                                                             pow(haplotypes[index + j].getprior(i)[n], sim))*/ ((1-sim) * haplotypes[index + j].posteriorwo[i][n] + sim * (neutraluncertain ? 0.5 : haplotypes[index + j].getprior(i)[n]));
                        datatarget[now][k + n] += datatarget[!now][k] * 0.5;
                    }
                }
            }
        }
        else
        {
            data[now].fill(1.f);
            datatarget[now].fill(1.f);
        }

        array<ratiotype, ploidy + 1> genotypebiasnow, target;
        double targetsum = 0;
        {
            //ratiotype sum = 0;
            genotypebiasnow.fill(0.f);
            target.fill(0.f);
            ratiotype sum = 0;
            for (int m = 0; m <= ploidy; m++)
            {
                ratiotype base = data[now][m];
                ratiotype targetbase = data[now][m];
                int counts[2] = {ploidy - m, m};
                // TODO WEAKENEPS DROPPEDz
                //if (genotypes[i] != -1 && counts[1] != genotypes[i]) base *= pow(std::max(domarkeps ? ourmap.otherepses[i] : 0.0f, epsothergeno) * (weakeneps ? (std::min(priors[j], priors[!j])) * 2 : 1.0f), abs(counts[1] - genotypes[i]));
                if (genotypes[i] != -1 && counts[1] != genotypes[i]) base *= pow(std::max(domarkeps ? ourmap.otherepses[i] : 0.0f, epsothergeno), abs(counts[1] - genotypes[i]));
                base *= genotypebias[counts[1]];
                targetbase *= genotypebias[counts[1]];

                for (int k = 0; k < 2; k++)
                {
                    if (reads[k] && !counts[k])
                    {
                        base *= 0;
                        targetbase *= 0;
                        continue;
                    }
                }                    

                int readsum = reads[0] + reads[1];
                for (int k = 0; k < 2; k++)
                {
                    for (int j = 0; j < reads[k]; j++)
                    {
                        base /= ploidy * 0.5;
                        targetbase /= ploidy * 0.5;
                        // Only one side of symmetry
                        if (!k)
                        {
                            base *= readsum - j;
                            base /= j + 1;
                            targetbase *= readsum - j;
                            targetbase /= j + 1;
                        }
                    }
                }
                // //base *= (m == 1) ? std::max(1 - haplotypes[index + 0].sim[i], 1e-10f) : 1;                

                sum += base;
                targetsum += targetbase;
                genotypebiasnow[m] = base;
                target[m] = targetbase;
            }
            sum += 1e-30f;
            sum = 1 / sum;
            targetsum += 1e-30f;
            targetsum = 1 / targetsum;
            for (int m = 0; m <= ploidy; m++)
            {
                genotypebiasnow[m] *= sum;
                target[m] *= targetsum;
            }
        }

        for (int m = 0; m < ploidy; m++)
        {
            double ourposteriormix = simposteriormix ? fabs(antisimposterior - haplotypes[index + m].sim[i]) :
                                    ((postmixred ? fabs(antisimposterior - haplotypes[index + m].sim[i]) : 1) * posteriormix);
            auto& priors = haplotypes[index + m].getprior(i);                                        

            if (newpostmix) ourposteriormix = 1 + haplotypes[index + m].sim[i] * (antiredcert ? certterm + priors[0] * priors[1] * certfactor : 1) * (-1 + posteriormix * (simredcert ? certterm + priors[0] * priors[1] * certfactor : 1));
            if (unknownred && reads[0] + reads[1] == 0 && genotypes[i] == -1) ourposteriormix = 0;

            if (redcertmix)
            {
                ourposteriormix *= certterm + priors[0] * priors[1] * certfactor;
            }
            ourposteriormix = 1;

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

            now = true;
            data[now].fill(0.f);
            data[now][0] = 1.0f;
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

            std::array<std::array<double, 2>, ploidy + 1> sums;
            for (int i = 0; i <= ploidy + 1; i++)
            {
                sums[i].fill(0.f);
            }
            for (int j = 0; j < 2; j++)
            {
                for (int a = 0; a < ploidy; a++)
                {
                    double base = data[now][a];
                    int counts[2] = {ploidy - 1 - a, a};
                    counts[j]++;

                    if ((burnin && mulpriorself) || (!burnin && propriorself)) base *= priors[j];
                    if (!burnin && antipriorself) base *= priors[!j];
                    if (!burnin && earlypost) base *= (selfposteriorwo && !nonwoearly) ?
                        (mixposteriors[m][j]) : haplotypes[index + m].posterior[i][j];
                    if ((burninassgn && burnin) || (postassgn && !burnin))
                    {
                        double notme = allnotme ? pow((ploidy - 1.0) / ploidy, reads[0] + reads[1]) :
                                                    (reads[j] ? pow((counts[j] - 1.0) / counts[j], reads[j]) : 1.0);
                        sums[counts[1]][j] += base * (1.0 - notme);
                        sums[counts[1]][0] += base * notme * ((earlypost && !burnin) ? ((selfposteriorwo && !selfpriorunass) ? haplotypes[index + m].posteriorwo[i][0] : haplotypes[index + m].posterior[i][0]) : 0.5);
                        sums[counts[1]][1] += base * notme * ((earlypost && !burnin) ? ((selfposteriorwo && !selfpriorunass) ? haplotypes[index + m].posteriorwo[i][1] : haplotypes[index + m].posterior[i][1]) : 0.5);
                    }
                    else
                    {
                        sums[counts[1]][j] += base;
                    }                        
                } 
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

            ratio[m] = 0;
            for (int i = 0; i <= ploidy; i++)
            {
                ratio[m] += sums[i][0] / (sums[i][0] + sums[i][1] + 1e-30f) * genotypebiasnow[i];
            }
        }
        ratiotype renormproportion = 0;
        for (int i = 0; i <= ploidy; i++)
        {
            renormproportion += genotypebiasnow[i] * (ploidy - i);
        }
        renormproportion /= ploidy;
        float speed = 0;
        
        if (!burnin)
            for (int i = 0; i <= ploidy; i++) { speed += pow(sqrt(genotypebiasnow[i]) - sqrt(target[i] * targetsum), 2); }
        else
            speed = 1;            
        //speed = sqrt(speed / 2);
        //for (int i = 0; i <= ploidy; i++) { speed += sqrt(genotypebiasnow[i]) * target[i] * targetsum); }
        updatenewpriors(i, ratio, speed /*renormproportion*/);
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
#pragma omp taskwait
    double extremesim = 0;
    int extrememarker = -1;
    for (int i = 0; i < genotypes.size(); i++)
    {
        if (!((genotypes[i] > 0 && genotypes[i] < ploidy) || (reads[i][0] && reads[i][1]))) continue;

        bool ok = true;
        for (int j = 0; j < ploidy; j++)
        {
            auto prior = haplotypes[index + j].getnewprior(i);
            if (prior[0] < 0.25 || prior[1] < 0.25) ok = false;
        }
        
        if (!ok) continue;
        
        double minsim = 1;
        for (int j = 0; j < ploidy; j++)
            for (int k = j + 1; k < ploidy; k++)
                if (haplotypes[index + j].crosssim[i][k] < minsim) minsim = haplotypes[index + j].crosssim[i][k];
        
        if (minsim > extremesim)
        {
            extrememarker = i;
            extremesim = minsim;
        }
    }

    if (extrememarker != -1)
    {
        int extremes[2] = {-1, -1};
        float extremevals[2] = {0, 0};
        for (int j = 0; j < ploidy; j++)
        {
            auto prior = haplotypes[index + j].getnewprior(extrememarker);
            for (int k = 0; k < 2; k++)
            {
                if (prior[k] > extremevals[k])
                {
                    extremevals[k] = prior[k];
                    extremes[k] = j;
                }
            }            
        }
        for (int k = 0; k < 2; k++)
        {
            auto& prior = haplotypes[index + extremes[k]].getnewprior(extrememarker);
            prior[k] = 1 - updeps;
            prior[!k] = updeps;
        }
        printf("Fixing marker %d at base index %d with similarity %f\n", extrememarker, index, extremesim);
    }
}

void zeroclasses()
{
    for (int i = 0; i < haplotypes.size(); i++)
    {
        haplotypes[i].classes.resize(haplotypes.size(), 0);
        for (int majorclass = 0; majorclass < nummajorclasses; majorclass++)
        {
            for (int j = 0; j < numclasses; j++)
            {
                haplotypes[i].initclassweights[majorclass][j] = (j == 0);
                for (int k = 0; k < numclasses; k++)
                {
                    haplotypes[i].classweights[majorclass][j][k] = ((j == 0) && (k == 0)) ? 1 : 0;
                }
            }
        }
    }
}

void normalizeclasses()
{
    for (int i = 0; i < haplotypes.size(); i++)
    {
        array<int, numclasses> counts{0};
        for (int& classval : haplotypes[i].classes)
        {
            counts[classval]++;
        }

        for (int majorclass = 0; majorclass < nummajorclasses; majorclass++)
        {
            for (int k = 0; k < numclasses; k++)
            {
                float sum = 0;
                for (int j = 0; j < numclasses; j++)
                {
                    if (counts[j]) sum += haplotypes[i].classweights[majorclass][j][k];
                }
                if (!sum) continue;
                sum = 1 / sum;
                for (int j = 0; j < numclasses; j++)
                {
                    haplotypes[i].classweights[majorclass][j][k] *= sum;
                }            
            }

            for (int j = 0; j < numclasses; j++)
            {
                float factor = 1;
                if (counts[j]) factor = 1.0f / counts[j];
                for (int k = 0; k < numclasses; k++)
                {                
                    haplotypes[i].classweights[majorclass][j][k] *= factor;
                }
                haplotypes[i].initclassweights[majorclass][j] *= factor;
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
        ind.maxshared.resize(ourmap.chromposes.size());
        ind.maxsharedid.resize(ourmap.chromposes.size(), 0);
        ind.samplehaplotypes(hapnum);
        hapnum += ploidy;
    }

    for (int h = basehaps; h < haplotypes.size(); h++)
    {
        auto& hap = haplotypes[h];
        for (int majorclass = 0; majorclass < nummajorclasses; majorclass++)
        {
            for (int fw = 0; fw < 2; fw++)
            {
                hap.renorm[majorclass][fw].resize(ourmap.chromposes.size());
            }
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
    globallelebias.resize(priors.size());
    #pragma omp taskloop
    for (int i = 0; i < priors.size(); i++)
    {
        genprob alleles = {0};
        int index = 0;
        for (genprob& g : priors[i])
        {
            alleles[0] += g[0];
            alleles[1] += g[1];
        }
        bool now = true;
        double sum = 0;        
        sum = alleles[0] + alleles[1];
        alleles[0] /= sum;
        alleles[1] /= sum;

        alleles[0] += updeps;
        alleles[1] += updeps;
        
        sum = 0;
        for (int i = 0; i < 2; i++)
        {
            sum += alleles[i];
        }
        sum = 1 / (sum + 1e-30f);
        for (int i = 0; i < 2; i++)            
        {
            alleles[i] *= sum;
        }
        sum = 0;
        for (int i = 0; i < 2; i++)
        {
            alleles[i] += updeps;
            alleles[i] = 1 / alleles[i];
            sum += alleles[i];
        }
        sum = 1 / sum;
        for (int j = 0; j < 2; j++)
        {
            globallelebias[i][j] = alleles[j] * sum;
        }
    }
    #pragma omp taskwait
    std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor + oneflip + relevel> fwbw[ploidy][nummajorclasses];
    int outer_tasks = omp_get_max_threads() / ploidy / 4 + 1;
    #pragma omp taskloop num_tasks(outer_tasks), private(hapnum, fwbw)
    //#pragma omp parallel for /*num_threads(16),*/ private(hapnum, fwbw)
    for (int i = 0; i < inds.size(); i++)
    {
        //#pragma omp parallel for num_threads(ploidy * 2), collapse(2), private(hapnum)
        //std::array<ArrayXXf, 2 + fullwo + 1 + nonsimfactor>* fwbw = ::fwbw;

        hapnum = basehaps + i * ploidy;
        for (int k = 0; k < ploidy; k++)
        {
            for (int majorclass = 0; majorclass < haplotypes[hapnum + k].herenummajor; majorclass++)
            {
                haplotypes[hapnum + k].fwbw[majorclass] = &fwbw[k][majorclass][0];
            }
        }
        for (int k = 0; k < ploidy; k++)
        {
            for (int majorclass = 0; majorclass < haplotypes[hapnum + k].herenummajor; majorclass++)
            {
                for (int fw = 0; fw < 2; fw++)
                #pragma omp task firstprivate(i, k, fw, hapnum)
                { 
                    individ& ind = inds[i];
                    haplotypes[hapnum + k].fwbw[majorclass][fw].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                    if (fw && fullwo)
                        haplotypes[hapnum + k].fwbw[majorclass][2].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                    if (fw /*&& ibdfactors*/)
                        haplotypes[hapnum + k].fwbw[majorclass][2 + fullwo].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                    if (fw /*&& ibdfactors*/ && nonsimfactor)
                        haplotypes[hapnum + k].fwbw[majorclass][2 + fullwo + nonsimfactor].resize(haplotypes.size() * 2, ourmap.chromposes.size());    
                    if (fw /*&& ibdfactors*/ && oneflip)
                        haplotypes[hapnum + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip].resize(haplotypes.size() * 2, ourmap.chromposes.size());
                    if (fw /*&& ibdfactors*/ && relevel && (k == 0 || multirelev))
                        haplotypes[hapnum + k].fwbw[majorclass][2 + fullwo + nonsimfactor + oneflip + relevel].resize(haplotypes.size() * 2, ourmap.chromposes.size());  
                    if (!burnin) haplotypes[hapnum + k].dofwbw(fw, ourmap, majorclass);
                }
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
            for (int majorclass = 0; majorclass < nummajorclasses; majorclass++)
            {
                haplotypes[hapnum + k].fwbw[majorclass] = nullptr;
            }
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
            //val = std::min(1.0f, val / (ploidy - 1) * (1 / (1 - val)));
            val = std::min(1.0f, val * (1 / (1 - val)));
        }
        if (val < 0) val = 0;
        ourmap.otherepses[i] = val; 
    }
}

double origstepsize = stepsize;

void doiter(int iter)
{
    if (iter == burniniters)
    {
        burnin = false;
        stepsize = origstepsize;
    }
    if (iter == setpostiter)
    {
        updallpriors = true;
        for (int i = 0; i < haplotypes.size(); i++)
        {
            for (int m = 0; m < haplotypes[i].posterior.size(); m++)
            {
                genprob& posts = haplotypes[i].posterior[m];
                float prisum = 0;
                genprob& priors = haplotypes[i].getprior(m);
                genprob& newpriors = haplotypes[i].getnewprior(m);
                for (float v : priors)
                    prisum += v;
                
                if (prisum)
                {
                    printf("prior already present for %d:%d\n", i, m);
                    continue;
                }

                float postsum = 0;
                for (float v : posts)
                    postsum += v;
                if (!postsum)
                {
                    printf("posterior not present for %d:%d\n", i, m);
                    continue;
                }

                for (int k = 0; k < priors.size(); k++)
                {
                    priors[k] = posts[k];
                    newpriors[k] = posts[k];
                }
                haplotypes[i].getanyprior(m) = 0.5;
                haplotypes[i].getnewanyprior(m) = 0.5;
            }
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
    likelihood = 0;
    doit();
}

#ifndef IMPUTATO_SKIP_MAIN
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
    for (auto& ind : inds)
    {
        if (ind.allowedrefs[0] == -1) continue;

        for (int i = 0; i < ind.genotypes.size(); i++)
        {
            if (ind.genotypes[i] == -1) continue;            

            int mins[2] = {2, 2};
            int maxs[2] = {0, 0};
            for (int j = 0; j < 2 * ploidy; j++)
            {
                int index = j / ploidy;
                int val = haplotypes[ind.allowedrefs[j]].getprior(i)[1] > 0.5;
                mins[index] -= 1 - val; //= std::min(val * 2, mins[index]);
                maxs[index] += val; //= std::max(val * 2, maxs[index]);
                mins[index] = std::max(0, mins[index]);
                maxs[index] = std::min(ploidy / 2, maxs[index]);
            }
            mins[0] += mins[1];
            maxs[0] += maxs[1];

            if (ind.genotypes[i] < mins[0] || ind.genotypes[i] > maxs[0])
            {
                printf("NOT WORKING %d:%d\t%d\t%lf\t(%d, %d)", indi, i,ind.genotypes[i], ourmap.otherepses[i], mins[0], maxs[0]);
                if (clearnonmendel) ind.genotypes[i] = -1;
                for (int j = 0; j < 2 * ploidy; j++)
                {
                    printf(" %d", haplotypes[ind.allowedrefs[j]].getprior(i)[1] > 0.5);
                }
                printf("\n");
            }
        }
        indi++;
    }
    initinds();
    zeroclasses();
    normalizeclasses();
    //inds.resize(2);
    burnin = true;
    stepsize = 0.2;
    for (int iter = 0; iter < itercount; iter++)
    {
        for (int i = 12; i < inds.size(); i += inds.size() - 1)
        {
            for (int j = 0; j < 200; j++)
            {
                if (!haplotypes[basehaps + i * ploidy].getanyprior(j) && !haplotypes[basehaps + i * ploidy].getanyprior(j + 1)) continue;

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
                printf("\t%.2f %d", inds[i].maxshared[j], inds[i].maxsharedid[j]);
                printf("\n");
            }
        }
        printf("Test! %d %lf\n", iter, likelihood);
        doiter(iter);
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
#endif