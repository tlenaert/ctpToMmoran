//
//  ctpgame.cpp
//  signalling
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include "ctpgame.hpp"
#include <cmath>
#include <iomanip>
#define FIXED_FLOAT(x) std::fixed <<std::setprecision(2)<<(x)


CtpGame::CtpGame(const CtpGame& other){
    _p1=other.p1();
    _p2=other.p2();
    _factor=other.factor();
    _length=other.length();
}

CtpGame& CtpGame::operator=(const CtpGame& other){
    _p1=other.p1();
    _p2=other.p2();
    _factor=other.factor();
    _length=other.length();
    return *this;
}


std::ostream&CtpGame::display(std::ostream& os) const {
    os << "("<< _length <<")[p1="  << FIXED_FLOAT(_p1) <<" ; p2="<< FIXED_FLOAT(_p2) <<" ; factor="<< FIXED_FLOAT(_factor) << "]";
    return os;
}

//Payoff as defined in the Rand paper
double CtpGame::payoff(unsigned round,bool which){
    //which is true = first mover, false = second mover
    if(round >=0 and round <= _length){
        double tmp1 = _p1* pow(_factor,round);
        double tmp2 = _p2* pow(_factor,round);;
        if(which)
            if(round%2==0)
                return tmp1;
            else return tmp2;
        else {
            if(round%2==0)
                return tmp2;
            else return tmp1;
        }
    }
    return NAN;
}

double CtpGame::payoff(unsigned rowA, unsigned rowB, unsigned colA, unsigned colB){//A and B define roles
    unsigned minimumA =min(rowA, colB); //row moves first
    double tmp1 = payoff(minimumA, true); // so payoff for row as first player
    unsigned minimumB =min(rowB, colA); // col moves first
    double tmp2 = payoff(minimumB, false); // so payoff for row as second player
    
    //need to adjust payoffs depending on who took the money when
    if(rowA != colB){
        if(rowA< colB && (rowA%2)!=0)
            tmp1 = payoff(minimumA+1, true);
        else if(colB < rowA && (colB%2)==0)
            tmp1 = payoff(minimumA+1, true);
    }
    if(rowB != colA){
        if(colA< rowB && (colA%2)!=0)
            tmp2 = payoff(minimumB+1, false);
        else if(rowB < colA && (rowB%2)==0)
            tmp2 = payoff(minimumB+1, false);
    }
    return (tmp1 + tmp2)/2.0; // return average payoff in both role A and role B
}



