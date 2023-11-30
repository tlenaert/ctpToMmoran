//
//  moran.cpp
//  signal
//
//  Created by Tom Lenaerts on 29/03/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include "moran.hpp"
#include <gsl/gsl_sf_exp.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "state.hpp"
#include <unordered_map>

#define FIXED_FLOAT(x) std::fixed <<setprecision(3)<<(x)

bool Moran::selectWithReplacement(Population* pop, RanGen* ran, unsigned amount, vector<unsigned>& result){
    for(unsigned i=0 ; i < amount ; i++){//need to check if the selected is present _distribution[selected]>0
        unsigned loc =ran->ranval(1, pop->pSize());
        unsigned selected = pop->selectNth(loc);
        result.push_back(selected);
    }
    return true;
}

bool Moran::selectWithoutReplacement(Population* pop, RanGen* ran, unsigned amount, vector<unsigned>& result){
    for(unsigned i=0 ; i < amount ; i++){ //need to check if the selected is present _distribution[selected]>0
        unsigned loc =ran->ranval(1, pop->pSize());
        while(std::find(result.begin(), result.end(), loc) != result.end())
            loc =ran->ranval(1, pop->pSize());
        unsigned selected = pop->selectNth(loc);
        if(pop->notEmpty(selected))
            result.push_back(selected);
        else {
            cout << "Selected empty strategy " << endl;
            exit(-1);
        }
    }
    return true;
}


double Moran::calcPayoff(Strategy* A,Strategy* B, RanGen* ran, double epsilon, unsigned repeats){ //identical to analytical code
    long double stochpayoff(0);
    for(unsigned i=0; i<repeats; i++){
        A->stochasticInferDecision(&_game, epsilon, ran);
        B->stochasticInferDecision(&_game, epsilon, ran);
//        cout << *first << "\t" << *second << "\t";
        long double tmp = _game.payoff(A->decisionR1(), A->decisionR2(), B->decisionR1(), B->decisionR2());
//        cout << tmp << endl;
        stochpayoff += tmp;
    }
    stochpayoff /=  double(repeats);
    return stochpayoff;
}



bool Moran::preCalcPayoff(StrategySpace* strats, RanGen* ran, double epsilon, unsigned repeats){
    if(_payoffs != NULL){
        gsl_matrix_free(_payoffs);
       
    }
    _payoffs = gsl_matrix_alloc(strats->size(), strats->size());
    gsl_matrix_set_zero(_payoffs);
    for(unsigned i=0 ; i< strats->size() ; i++){
        Strategy* A = (*strats)[i];
        for(unsigned j=0 ; j< strats->size() ; j++){
            Strategy* B = (*strats)[j];
            double tmp = calcPayoff(A, B, ran, epsilon, repeats);
            gsl_matrix_set(_payoffs, i, j, tmp);
        }
    }
    return true;
}

void Moran::printPayoffs(StrategySpace* space){
    cout << "\t\t\t";
    for (unsigned i = 0; i < _payoffs->size1; i++) {
        cout << *((*space)[i]) << "\t";
    }
    cout << endl;
    for (unsigned i = 0; i < _payoffs->size1; i++) {
        cout << *((*space)[i]) << "\t";
        for (size_t j = 0; j < _payoffs->size2; j++) {
            cout << FIXED_FLOAT(gsl_matrix_get(_payoffs, i, j));
            if( j < (_payoffs->size2 - 1))
                cout<<"\t";
        }
        cout << endl;
    }
}


double Moran::playAgainstall(unsigned other, Population* pop){
    double fitness = 0.0;
    double total = 0;
    for(unsigned i=0 ; i< pop->numStrats(); i++){
        double numb =  (double)(*pop)[i];
        total += numb;
        if(i != other){
            fitness += numb * gsl_matrix_get(_payoffs, other, i);
        }
        else {
            fitness += (numb-1.0)* gsl_matrix_get(_payoffs, other, i);
        }
    }
    return (fitness / double(total-1));
}

double Moran::playAgainstall(Strategy* other, Population* pop, StrategySpace* space, RanGen* ran,double epsilon, unsigned repeats){
    double fitness = 0.0;
    double total = 0;
    for(unsigned i=0 ; i< pop->numStrats(); i++){
        Strategy* elm = (*space)[i];
        double numb =  (double)(*pop)[i];
        total += numb;
        if(*elm != *other){
            fitness += numb * calcPayoff(other, elm, ran, epsilon, repeats);
        }
        else {
            fitness += (numb-1.0)* calcPayoff(other, elm, ran, epsilon, repeats);
        }
    }
    return (fitness / double(total-1));
}

bool Moran::imitate(RanGen *ran, double first, double second){ // fermifunction
    double result = 1.0 / (1.0 + gsl_sf_exp((_beta) *(first-second)));
    return (ran->randouble() < result);
}

bool Moran::moranStepPairwise(unsigned& step, Population* pop, RanGen* ran, StrategySpace* space, unordered_map<unsigned long, State>& storage, double epsilon, unsigned repeats){
    vector<unsigned> selected;
    unsigned toselect=2;
    selectWithoutReplacement(pop, ran, toselect, selected);
    
    vector<double> fitness;
    for (unsigned i=0; i < toselect; i++){
//        double value = playAgainstall((*space)[selected[i]], pop, space, events,ran, epsilon, repeats);
        double value = playAgainstall(selected[i], pop); // when fitness is precalculated.
        fitness.push_back(value);
    }
    
    if(imitate(ran,fitness[0], fitness[1])){
        if(ran->randouble() < _mut){
            unsigned loc = ran->ranval(0,pop->numStrats()-1);
            pop->swap(selected[0],loc,1);
        }
        else {
            pop->swap(selected[0],selected[1],1);
        }
    }
    else {// 1 imitates 0
        if(ran->randouble() < _mut){
            unsigned loc = ran->ranval(0,pop->numStrats()-1);
            pop->swap(selected[1],loc,1);
        }
        else {
            pop->swap(selected[1],selected[0],1);
        }
    }

    //collect population state information
    State tmp(pop, space);
    auto located = storage.find(tmp.index());
    if(located == storage.end()){
        located = storage.insert(make_pair(tmp.index(), tmp)).first;
    }
    else {
        located->second.increase(1);
    }
    
    int test = pop->converged();
    if(_mut> 0 && test != -1){ //jump froward to the next mutation once converged
//        cout << *pop << "\t";
        unsigned jump = ran->rangeometric(_mut);
	if(jump < (_iterations - step)){
        	// force mutation to new state
        	unsigned loc = ran->ranval(0,pop->numStrats()-1);
//        	cout << *pop << endl;
        	while(loc == test)
            		loc = ran->ranval(0,pop->numStrats()-1);
        	pop->swap(test,loc,1);
//        	cout << *pop << endl;
        	//update state
        	located->second.increase(jump-1);
//        	cout << located->second << endl;
        	step += jump;
//        	cout << "jumped " << (jump-1) << endl;
		State tmpp(pop, space);
            	auto located = storage.find(tmpp.index()); //change this to a string key; concatenation of all values in distribution.
            	if(located == storage.end()){
                	located = storage.insert(make_pair(tmpp.index(), tmpp)).first;
            	}
            	else {
                	located->second.increase(1);
            	}
	}
        else {
          	located->second.increase(_iterations-step-1);
		step=_iterations;
		cout << "Wanted to jump beyond number of iterations; " << (step+jump) << endl;
	}
    }
    selected.clear();
    return true;
}


void Moran::addStratDistribution(gsl_vector* total){
    gsl_vector_add(total, _composition);
}


//bool compareState(pair<unsigned long, State>& a, pair<unsigned long, State>& b){
//    return a.second.visits() < b.second.visits();
//}

struct compareStates{
    bool operator()(pair<unsigned long, State>& a, pair<unsigned long, State>& b){
        return a.second.visits() > b.second.visits();
    }
};


bool Moran::execute(Population* pop, RanGen* ran, StrategySpace* space, unsigned iterations, double epsilon, unsigned repeats){
    _psize=pop->pSize();
    _iterations=iterations;
    
    //states
    unordered_map<unsigned long, State> storage;

    //launch moran process.
    unsigned steps=0;
    for(; steps < iterations ; steps++){
        moranStepPairwise(steps, pop, ran, space, storage, epsilon, repeats);
    }

    //initialise strategies vector
    if(_composition !=NULL)
        gsl_vector_free(_composition);
    _composition = gsl_vector_alloc(space->size());
    gsl_vector_set_zero(_composition);


    // average appearance of startegy over all states is stored in _composition
    double total = 0.0;
    unsigned visits=0;
    for(auto iter = storage.begin(); iter != storage.end(); iter++ ){
        State tmp = iter->second;
        visits+= tmp.visits();
        for(unsigned i = 0; i < tmp.size(); i++){
            if(tmp[i]>0){
                double val = double(tmp.visits()) * (double(tmp[i])/double(_psize));
                gsl_vector_set(_composition, i, gsl_vector_get(_composition, i)+val);
            }
        }
        total+=1.0;
    }
    cout << "states collected = " << round(total) << ", total visits " << visits << ", steps in moran " << round(steps) << endl;
    double test=0.0;
    
    for(unsigned i=0; i < _composition->size; i++){
        double val = gsl_vector_get(_composition, i);
        gsl_vector_set(_composition, i, (val / (double(steps))));
        test += (val / (double(steps)));
    }
    cout << "Sum distrib = " << test << endl;
    storage.clear();
    return true;
}
