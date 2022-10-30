//
//  ctpgame.hpp
//  signalling
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#ifndef ctpgame_hpp
#define ctpgame_hpp

#include <string.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

class CtpGame {
public:
    CtpGame(double p1, double p2, double factor, unsigned length):_p1(p1), _p2(p2), _factor(factor),_length(length){};
    CtpGame(const CtpGame& other);
    
    double p1() const {return _p1;}
    bool setP1(double val) {_p1=val;return (_p1==val);}
    double p2() const {return _p2;}
    bool setP2(double val) {_p2=val;return (_p2==val);}
    double factor() const {return _factor;}
    bool setFactor(double val) {_factor=val;return (_factor==val);}
    unsigned length() const {return _length;}
    bool setLength(unsigned val) {_length=val;return (_length==val);}
    
    virtual double payoff(unsigned rounds, bool which);
    virtual double payoff(unsigned rowA, unsigned rowB, unsigned colA, unsigned colB);

    friend ostream & operator<<(ostream &o, CtpGame& g){return g.display(o);}
    virtual CtpGame& operator=(const CtpGame& other);

protected:
    virtual ostream& display(ostream& os) const ;

    
    bool check();
    double _p1;
    double _p2;
    double _factor;
    unsigned _length;
    
    //4 game data mc
};


#endif /* ctpgame_hpp */
