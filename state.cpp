//
//  state.cpp
//  norms
//
//  Created by Tom Lenaerts on 30/04/2021.
//

#include "state.hpp"
#include <algorithm>
#include <sstream>


unsigned long State::binomialCoeff(unsigned n, unsigned k) {
    unsigned long res = 1;
    // Since C(n, k) = C(n, n-k)
    if (k > n - k) k = n - k;
    // Calculate value of [n * (n-1) * ... * (n-k+1)] / [k * (k-1) * ... * 1]
    for (size_t i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}


unsigned long State::starsBars(unsigned stars, unsigned bins) {
    return binomialCoeff(stars + bins - 1, stars);
}

void State::calcIndex(unsigned psize) {
    //from EGTtools
    _index = 0;
    unsigned remaining = psize;
    for (unsigned i = 0; i < _distribution.size() - 1; ++i) {
        unsigned h = remaining;
        while (h > _distribution[i]) {
            _index += starsBars(remaining - h, _distribution.size() - i - 1);
            --h;
        }
        if (remaining == _distribution[i])
            break;
        remaining -= _distribution[i];
    }
}

void State::indexToState(vector<unsigned>& reconstruct) {
    //from EGTtools
    auto remaining = _popsize;
    unsigned nb_strategies = (unsigned)_distribution.size();
    reconstruct.clear();reconstruct.resize(nb_strategies);
    if(_index > 0){
        unsigned long i = _index;
        for (unsigned a = 0; a < nb_strategies; ++a) {
            // reset the state container
            reconstruct[a] = 0;
            for (unsigned j = remaining; j > 0; --j) {
                unsigned long count= starsBars(remaining - j, nb_strategies - a - 1);
                if (i >= count) {
                    i -= count;
                } else {
                    reconstruct[a] = j;
                    remaining -= j;
                    break;
                }
            }
        }
    }
}

State::State(Population* pop, StrategySpace* space){
    _distribution.clear();
    _popsize = pop->pSize();
    _visits = 1;
    for(unsigned i=0; i < pop->numStrats(); i++)
        _distribution.push_back((*pop)[i]);
    
    calcIndex(_popsize);
}
    

State::State(const State& other){
    _distribution.clear();
    _popsize = other.psize();
    for(unsigned i=0; i < other.size(); i++)
        _distribution.push_back(other[i]);
    _visits = other.visits();
    _index = other.index();
}

State& State::operator=(const State& other){
    _distribution.clear();
    _popsize = other.psize();
    for(unsigned i=0; i < other.size(); i++)
        _distribution.push_back(other[i]);
    _visits = other.visits();
    return *this;
}

bool State::compareID(){
    vector<unsigned> rtest;
    this->indexToState(rtest);

    if(_distribution.size() != rtest.size())
        return false;
    for(unsigned i=0; i < _distribution.size(); i++){
        if (_distribution[i] != rtest[i])
            return false;
    }
    return true;

}

bool State::operator==(const State& other) const{
    if(_distribution.size() != other.size() || _popsize != other.psize())
        return false;
    if(_index != other.index())
        return false;
    for(unsigned i=0; i < _distribution.size(); i++){
        if (_distribution[i] != other[i])
            return false;
    }
    return true;
}

bool State::operator!=(const State& other) const{
    return !(*this == other);
}


ostream& State::display(std::ostream& os) const {
    os << "[" << _index << " : " << _visits;
    os << "]";
    return os;
}

unsigned State::operator[](unsigned pos) const{
    if(pos >=0 && pos < _distribution.size())
        return _distribution[pos];
    return (unsigned)NAN;
}
