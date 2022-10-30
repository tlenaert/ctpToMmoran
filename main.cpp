//
//  main.cpp
//  ctpbeliefmoran
//
//  Created by Tom Lenaerts on 14/05/2022.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <chrono>
#include <gsl/gsl_matrix.h>

#include "ctpgame.hpp"
#include "ctpdata.hpp"
#include "rangen.h"
#include "strategy.hpp"
#include "population.hpp"
#include "moran.hpp"


using namespace std;
using namespace std::chrono;
#define FIXED_FLOAT(x) std::fixed <<setprecision(9)<<(x)


void extractBeliefDistrubution(StrategySpace* strategies, gsl_vector* distribution, ofstream& of, double runs, double beta, double eps){
    map<string,double> results;
    for(unsigned i=0; i < strategies->size();i++){
        double tmp = (gsl_vector_get(distribution, i)/runs);
        Strategy* elm = (*strategies)[i];
        stringstream ss;
        ss << "(" << elm->beliefR1() << "," << elm->beliefR2() <<")";
        map<string,double>::iterator found =results.find(ss.str());
        if(found != results.end()){
            found->second +=tmp;
        }
        else results[ss.str()] = tmp;
    }

    map<string,double>::iterator start = results.begin();
    map<string,double>::iterator stop = results.end(); // I assume the iteration is always in the same order
    
    of << beta << "\t" << eps <<"\t";
    unsigned numb= (unsigned)results.size();
    unsigned count=0;
    while(start != stop){
        cout << start->first << "\t";
        of << start->second;
        if(count < (numb-1))
            of << "\t";
        start++;
        count++;
    }
    of << endl;
    cout << endl;
}

void extractLevelDistrubution(StrategySpace* strategies, gsl_vector* distribution, ofstream& of, double runs, unsigned levels, double beta, double eps){
    map<unsigned,vector<double>> results;

    for(unsigned i=0; i < strategies->size();i++){
        double tmp = (gsl_vector_get(distribution, i)/runs);
        Strategy* elm = (*strategies)[i];
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        if(found != results.end()){
            vector<double> data = found->second;
            data[elm->beliefR1()] += tmp;  //symmetric beliefs
            found->second = data;
        }
        else {
            vector<double> data(levels, 0.0);
            data[elm->beliefR1()] += tmp;
            results[elm->level()] = data;
        }
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();

    while (start!=stop){
        unsigned lev = start->first;
        of << beta << "\t" << eps << "\t" << lev << "\t";
        double total = 0;
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            of << data[i];
            total += data[i];
            if (i < (data.size()-1))
                of << "\t";
        }
        of << "\t" << total << endl;
        start++;
    }
    of << endl;
}

void extractStepsDistrubution(StrategySpace* strategies, gsl_vector* distribution, ofstream& of, double runs, RanGen* ran, CtpGame* game, double epsilon, unsigned repeats, double beta, unsigned levels){
    vector<double> results(levels,0.0);
    for(unsigned i=0; i < strategies->size();i++){
        double tmp = (gsl_vector_get(distribution, i)/runs);
        Strategy* elm = (*strategies)[i];
        for(unsigned j=0; j < repeats; j++){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(game, epsilon, ran);
            second.stochasticInferDecision(game, epsilon, ran);
            unsigned whenp1=first.decisionR1();
            if(whenp1%2!=0 && whenp1 != game->length())
                whenp1+=1;
            unsigned whenp2=second.decisionR2();
            if(whenp2%2==0 && whenp2 != game->length())
                whenp2+=1;
            if(whenp1 <= whenp2)
                results[whenp1] += tmp;
            else results[whenp2] += tmp;
        }
    }
    of << beta << "\t" << epsilon << "\t";
    double sum=0;
    for (unsigned j=0; j < results.size(); j++){
        of << (results[j]/double(repeats));
        sum+=results[j];
        if (j < (results.size()-1))
            of << "\t";
    }
    of << endl;
}

void extractDecisionDistrubution(StrategySpace* strategies, gsl_vector* distribution, ofstream& of, double runs, RanGen* ran, CtpGame* game, double epsilon, unsigned repeats, double beta, unsigned levels){
    map<unsigned,vector<double>> results;

    for(unsigned i=0; i < strategies->size();i++){
        double tmp = gsl_vector_get(distribution, i);
        Strategy* elm = (*strategies)[i];
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        vector<double> collect(levels,0.0);
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(game, epsilon, ran);
            second.stochasticInferDecision(game, epsilon, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            map<unsigned,double>::iterator found;
            unsigned index=whenp1;
            if(whenp1 > whenp2)
                index=whenp2;
            collect[index] += tmp;
        }
        if(found != results.end()){
            vector<double> data = found->second;
            for(unsigned iter=0; iter < levels; iter++){
                double val = collect[iter]/double(repeats);
                data[iter] += val;
            }
            found->second = data;
        }
        else {
            for(unsigned iter=0; iter < levels; iter++){
                collect[iter]/=double(repeats);
            }
            results[elm->level()] = collect;
        }
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();

    while (start!=stop){
        unsigned lev = start->first;
        of << beta << "\t" << epsilon << "\t" << lev << "\t";
        double total = 0;
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            of << data[i];
            total += data[i];
            if (i < (data.size()-1))
                of << "\t";
        }
        of << "\t" << total << endl;
        start++;
    }
    of << endl;
}

void extractMisbeliefDistrubution(StrategySpace* strategies, gsl_vector* distribution, ofstream& of, double runs, RanGen* ran, CtpGame* game, double epsilon, unsigned repeats, double beta, unsigned levels){
    map<unsigned,vector<double>> results;
    unsigned size = ((2*levels)-1);
    unsigned middle = (unsigned)floor(size/2);

    for(unsigned i=0; i < strategies->size();i++){
        Strategy* elm = (*strategies)[i];
        double tmp = gsl_vector_get(distribution, i);
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        vector<double> collect(size,0.0);
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(game, epsilon, ran);
            second.stochasticInferDecision(game, epsilon, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            map<unsigned,double>::iterator found;
            unsigned index=whenp1;
            if(whenp1 > whenp2)
                index=whenp2;
            int misbelief = (elm->beliefR1() - index); // symmetric belief
            collect[middle+misbelief] += tmp;
        }
        double test = 0;
        if(found != results.end()){
            vector<double> data = found->second;
            for(unsigned iter=0; iter < data.size(); iter++){
                double val = collect[iter]/double(repeats);
                data[iter] += val;
                test+=val;
            }
            found->second = data;
        }
        else {
            for(unsigned iter=0; iter < collect.size(); iter++){
                collect[iter]/=double(repeats);
                test+=collect[iter];
            }
            results[elm->level()] = collect;
        }
//        cout << *elm << "\t" << tmp << "\t" << test << endl;
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();
    while (start!=stop){
        unsigned lev = start->first;
        of << beta << "\t" << epsilon << "\t" << lev << "\t";
        double pos(0), correct (0), neg(0), scaled(0), sum(0);
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            sum+=data[i];
            if(i < middle) neg += data[i];
            if(i > middle) pos += data[i];
            if (i== middle) correct+=data[i];
            if(i < middle)
                scaled -= (data[i] * (levels-i-1));
            else scaled += (data[i] * (i-levels+1));
        }

        of << (correct/sum) << "\t" << (neg/sum) << "\t" << (pos/sum) << "\t" << scaled << endl;  //prints now the fraction not scaled to stationary distribution
        start++;
    }
    of << endl;
}

void runMoranSimulations(unsigned psize, double beta, double epsilon, double mut, unsigned levels, unsigned maxlevel, unsigned repeats, CtpGame& game, CtpData& data, unsigned runs, unsigned iterations, ofstream& bf, ofstream& lf, ofstream& sf, ofstream& df, ofstream& mf, RanGen& ran){

    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);

    Population pop;
    RandomSimplexPopulationFactory bob(strategies.size(), &ran); // test this works correctly

    Moran process(psize, beta, mut, game);

    process.preCalcPayoff(&strategies,&ran, epsilon, repeats);
    cout << strategies << endl;


    gsl_vector* stratdistrib = gsl_vector_alloc(strategies.size()); // no + all levels
    gsl_vector_set_zero(stratdistrib);


    for(unsigned iter=0; iter < runs; iter++){
        high_resolution_clock::time_point start = high_resolution_clock::now();
        cout<< "[NEW RUN : " << iter <<"]" << endl;
        if(bob.createPopulation(psize, &pop)){
            cout << pop <<  pop.sum() << endl;
            process.execute(&pop, &ran, &strategies, iterations, epsilon, repeats);
            process.addStratDistribution(stratdistrib);
        }
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by function: "
             << duration.count() << " s" << endl;
    }


    extractBeliefDistrubution(&strategies, stratdistrib, bf, runs, beta, epsilon);
    extractStepsDistrubution(&strategies, stratdistrib, sf, runs, &ran, &game, epsilon, repeats, beta,levels+1);
    extractLevelDistrubution(&strategies, stratdistrib, lf, runs, levels+1, beta, epsilon);
    extractDecisionDistrubution(&strategies, stratdistrib, df, runs, &ran, &game, epsilon, repeats, beta,levels+1);
    extractMisbeliefDistrubution(&strategies, stratdistrib, mf, runs, &ran, &game, epsilon, repeats, beta,levels+1);

    gsl_vector_free(stratdistrib);
}


int main(int argc, char * argv[]) {
    unsigned length = 4;
    double first = 0.4;
    double second = 0.1;
    double factor= 2.0;
    unsigned levels = length;
    unsigned maxlevel=4;
    double epsilon=0.18;
    unsigned repeats=50000;
    unsigned psize = 500;
    double betas=0.3;
    double mut = 0.0;
    unsigned runs=100;
    unsigned iterations = 10000000;
    double cost = 0.0;

    //output filenames
    stringstream ss1;
    ss1 << "./ctplevelsL" << maxlevel << ".txt";
    string lfname(ss1.str());
    stringstream ss2;
    ss2 << "./ctpbeliefsL" << maxlevel << ".txt";
    string bfname(ss2.str());
    stringstream ss3;
    ss3 << "./ctpstepsL" << maxlevel << ".txt";
    string sfname(ss3.str());
    stringstream ss4;
    ss4 << "./ctpdecisionsL" << maxlevel << ".txt";
    string dfname(ss4.str());
    stringstream ss5;
    ss5 << "./ctpmisbeliefL" << maxlevel << ".txt";
    string mfname(ss5.str());

   
    
    cout << "Settings: ";
    cout << "psize(" << psize << "),  Beta-s(" << betas << "), Eps(" << epsilon << "), ";
    cout << "Length(" << length <<"), Factor(" << factor <<"), First(" << first <<"), Second(" << second <<"), Levels(" << maxlevel <<"), " ;
    cout <<  "Mut(" << mut << "), repeats("<<repeats<<"), " ;
    cout << "Runs(" << runs << "),  Iterations(" << iterations << "), ";
    cout << "LFname ("<< lfname <<"), ";
    cout << "BFname ("<< bfname <<"), ";
    cout << "SFname ("<< sfname <<"), ";
    cout << "DFname ("<< dfname <<"), ";
    cout << "MFname ("<< mfname <<"), " << endl;

    // prepare the games, rangen, data and file variables to execute the Moran simulation.
    CtpGame game(first, second, factor, length);
    RanGen ran;
    
    CtpData data;
    CtpEntry MKP4avg(4, {0.071, 0.356, 0.370, 0.153, 0.049});
    data.add("MKP4avg", &MKP4avg);
    CtpEntry MKP6avg(6, {0.007, 0.064, 0.199, 0.384, 0.253, 0.078, 0.014});
    data.add("MKP6avg", &MKP6avg);
    cout << data;

    ofstream lf(lfname);
    ofstream sf(sfname);
    ofstream bf(bfname);
    ofstream df(dfname);
    ofstream mf(mfname);

    runMoranSimulations(psize, betas, epsilon, mut, levels, maxlevel, repeats, game, data, runs, iterations, bf, lf, sf, df, mf, ran);


    lf.close();
    sf.close();
    bf.close();
    df.close();
    mf.close();
    return 0;
}
