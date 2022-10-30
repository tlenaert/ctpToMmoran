//
//  population.cpp
//  states
//
//  Created by Tom Lenaerts on 30/11/2018.
//  Copyright © 2018 Tom Lenaerts. All rights reserved.
//

#include "population.hpp"
#include <algorithm>

bool Population::reset(unsigned strategies){
    if(_distribution!=NULL)
        delete[] _distribution;
    _distribution= new unsigned[strategies];
    _distrib_size = strategies;
    for(unsigned i=0; i < _distrib_size; i++)
        _distribution[i]=0;
    _psize=0;
    return true;
}

void Population::set(unsigned loc, unsigned amount) {
    if(loc >=0 && loc < _distrib_size)
        _distribution[loc]=amount;
}
void Population::increase(unsigned loc, unsigned amount) {
    if(loc >=0 && loc < _distrib_size && _distribution[loc] < _psize)
        _distribution[loc]+=amount;
}
void Population::decrease(unsigned loc, unsigned amount) {
    if(loc >=0 && loc < _distrib_size && _distribution[loc] > 0)
        _distribution[loc]-=amount;
}

void Population::swap(unsigned from, unsigned to, unsigned amount){
    if(from >=0 && from < _distrib_size && to >=0 && to < _distrib_size){
        _distribution[from]-= amount;
        _distribution[to]+= amount;
    }
}

bool Population::notEmpty(unsigned loc){
    return (loc >=0 && loc < _distrib_size && _distribution[loc]!=0);
}

unsigned Population::operator[](unsigned pos) const{
    if (pos >=0 && pos < _distrib_size)
        return _distribution[pos];
    return (unsigned)NAN;
}


unsigned Population::selectNth(unsigned amount){
    unsigned total = 0;
    int iter=-1;
    do{
        iter +=1;
        total += _distribution[iter];
    }
    while(total < amount);
    return iter;
}




int Population::converged(){
    for(unsigned i=0; i < _distrib_size; i++){
        if (_distribution[i] == pSize()){
            return i;
        }
    }
    return -1;
}

unsigned Population::sum(){
    unsigned total = 0;
    for(unsigned i=0; i < _distrib_size; i++){
        total += _distribution[i];
    }
    return total;
}


ostream& Population::display(std::ostream& os) const {
    os<<"{";
    for (unsigned iter=0 ; iter < _distrib_size; ++iter){
        os << iter << " : " << _distribution[iter] << endl;
    }
    os<<"}" << endl;
    return os;
}




bool UniformRandomPopulationFactory::createPopulation(unsigned psize, Population* pop){
    pop->reset(_strategies);
    for (unsigned i=0; i < psize; i++){
        unsigned selected = _ran->ranval(0,_strategies-1);
        pop->increase(selected,1);
    }
    pop->setPSize(psize);
        
    return true;
}

bool EqualSplitPopulationFactory::createPopulation(unsigned psize, Population* pop){
    pop->reset(_strategies);
    unsigned splits=floor(psize/ double(_strategies));
    unsigned remainder = psize % _strategies;
    if(remainder != 0){
        psize += (_strategies-remainder);
        splits=floor(psize/ double(_strategies));
    }
    for (unsigned i=0; i < _strategies; i++){
        pop->set(i, splits);
    }
    pop->setPSize(psize);

    return true;
}

bool RandomSimplexPopulationFactory::createPopulation(unsigned psize, Population* pop){
    //based on https://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf
    pop->reset(_strategies);
    unsigned remainder = psize % _strategies;
    if(remainder != 0){
        psize += (_strategies-remainder);
    }
    
    vector<unsigned> rannums;
    rannums.push_back(0);
    for (unsigned i=1; i < _strategies; i++){
        unsigned val = _ran->ranval(1, psize-1);
        while(find(rannums.begin(),rannums.end(), val) != rannums.end())
            val = _ran->ranval(1, psize-1);
        rannums.push_back(val);
    }
    rannums.push_back(psize);
    sort(rannums.begin(), rannums.end());
    unsigned check = 0;
    for (unsigned i=0; i < (rannums.size()-1); i++){
        check += (rannums[i+1]-rannums[i]);
        pop->set(i, rannums[i+1]-rannums[i]);
    }
    if(check != psize){
        cout << "Wrong population size, check is " << check << " and should be " << psize << endl;
        exit(-1);
    }
    pop->setPSize(psize);

    return true;
}
