//
//  moran.hpp
//  signal
//
//  Created by Tom Lenaerts on 29/03/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#ifndef moran_hpp
#define moran_hpp

#include <stdio.h>
#include "population.hpp"
#include "ctpgame.hpp"
#include "strategy.hpp"
#include "population.hpp"
#include <vector>
#include "rangen.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <unordered_map>
#include "state.hpp"

class Moran {
public:
    Moran(unsigned psize, double beta, double mut,  CtpGame& game):_psize(psize),  _beta(beta), _game(game), _mut(mut) {
        _payoffs = NULL;
        _composition=NULL;
        _iterations = 0;
    };
    
    ~Moran(){
        if(_composition != NULL)
            gsl_vector_free(_composition);
        if(_payoffs != NULL)
            gsl_matrix_free(_payoffs);
    }
    
    bool preCalcPayoff(StrategySpace* strats,RanGen* ran, double epsilon, unsigned repeats);
    void printPayoffs(StrategySpace* space);

    bool execute(Population*, RanGen*, StrategySpace*, unsigned, double, unsigned);
    void addStratDistribution(gsl_vector* total);

protected:
    bool selectWithReplacement(Population* pop, RanGen* ran, unsigned amount, vector<unsigned>& result);
    bool selectWithoutReplacement(Population* pop, RanGen* ran, unsigned amount, vector<unsigned>& result);
    double playAgainstall(unsigned other,Population* pop);
    double playAgainstall(Strategy* other,Population* pop, StrategySpace* space, RanGen* ran, double epsilon, unsigned repeats);
    bool imitate(RanGen *ran, double first, double second);
    bool moranStepPairwise(unsigned& step,Population* pop, RanGen* ran,StrategySpace* space,unordered_map<unsigned long, State>& storage, double epsilon, unsigned repeats);
    double calcPayoff(Strategy* A,Strategy* B, RanGen* ran, double eps, unsigned repeats);
    
    unsigned _psize;
    double _beta;
    CtpGame _game;
    double _mut;
    unsigned _iterations;
    
    gsl_vector* _composition;
    gsl_matrix* _payoffs;
};


#endif /* moran_hpp */
