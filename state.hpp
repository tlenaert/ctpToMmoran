//
//  state.hpp
//  norms
//
//  Created by Tom Lenaerts on 30/04/2021.
//

#ifndef state_hpp
#define state_hpp

#include <stdio.h>
#include <string.h>
#include "population.hpp"

using namespace std;

class State {
public :
    State(){
        _distribution.clear();
        _popsize = 0;
        _index = 0;
        _visits = 0;
    };
    State(Population* pop,StrategySpace* space);
    
    ~State(){};
    
    State(const State& other);
    State& operator=(const State& other);
    
    unsigned operator[](unsigned pos) const;

    unsigned long index() const {return _index;}
    void indexToState(vector<unsigned>& reconstruct);
    bool compareID();

    unsigned size() const {return (unsigned int)_distribution.size();}
    unsigned psize() const {return _popsize;}
    unsigned visits() const {return _visits;}
    void increase(unsigned val) {_visits +=val;}
    void reset() {_visits = 1;}

    bool operator==(const State& other) const;
    bool operator!=(const State& other) const;

    friend ostream & operator<<(ostream &o, State& s){return s.display(o);}
   
protected:
    virtual ostream& display(ostream& os) const ;
    void calcIndex(unsigned psize);
    unsigned long starsBars(unsigned stars, unsigned bins);
    unsigned long binomialCoeff(unsigned n, unsigned k);
    
    vector<unsigned> _distribution;
    unsigned long _index;
    unsigned _popsize;
    unsigned _visits;
};


#endif /* state_hpp */
