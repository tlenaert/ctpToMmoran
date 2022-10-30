//
//  population.hpp
//  states
//
//  Created by Tom Lenaerts on 30/11/2018.
//  Copyright © 2018 Tom Lenaerts. All rights reserved.
//

#ifndef population_hpp
#define population_hpp

#include <stdio.h>
#include <vector>

#include "rangen.h"
#include <iostream>
#include <vector>
#include "strategy.hpp"
#include <sstream>
#include <fstream>

using namespace std;

class Population {
public:
    Population(){
        _distribution = NULL;
    };

    ~Population(){
        if(_distribution!=NULL)
            delete[] _distribution;
    }
        
    bool reset(unsigned strategies);
    
    unsigned pSize() const {return _psize;}
    void setPSize(unsigned psize) {_psize=psize;}
    unsigned numStrats() const {return _distrib_size;}
    void setStratNum(unsigned dsize) {_distrib_size=dsize;}

    unsigned selectNth(unsigned amount);
    
    void set(unsigned loc, unsigned amount);
    void increase(unsigned loc, unsigned amount);
    void decrease(unsigned loc, unsigned amount);
    void swap(unsigned from, unsigned to, unsigned amount);
    bool notEmpty(unsigned loc);
    int converged();
    unsigned sum();
    
    unsigned operator[](unsigned pos) const;

    friend std::ostream & operator<<(std::ostream &o, Population& p){return p.display(o);}

protected:
    virtual std::ostream& display(std::ostream& os) const ;
    unsigned *_distribution;
    unsigned _distrib_size;
    unsigned _psize;
};


class PopulationFactory{
public:
    virtual bool createPopulation(unsigned, Population*)=0;
    
};

class UniformRandomPopulationFactory : public PopulationFactory{
public:
    UniformRandomPopulationFactory(unsigned strategies, RanGen* ran):_strategies(strategies), _ran(ran){};
    ~UniformRandomPopulationFactory(){
        _ran=NULL;
    }
    
    bool createPopulation(unsigned psize, Population* pop);
    
protected:
    unsigned _strategies;
    RanGen *_ran;
};


class EqualSplitPopulationFactory : public PopulationFactory{
public:
    EqualSplitPopulationFactory(unsigned strategies):_strategies(strategies){};
    bool createPopulation(unsigned psize, Population* pop);
    
protected:
    unsigned _strategies;
};


class RandomSimplexPopulationFactory : public PopulationFactory{
public:
    RandomSimplexPopulationFactory(unsigned strategies, RanGen* ran):_strategies(strategies), _ran(ran){};
    ~RandomSimplexPopulationFactory(){
        _ran=NULL;
    }
    
    bool createPopulation(unsigned psize, Population* pop);
    
protected:
    unsigned _strategies;
    RanGen *_ran;
};

#endif /* population_hpp */
