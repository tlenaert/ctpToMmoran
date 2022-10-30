//
//  ctpdata.hpp
//  ctpbeliefs
//
//  Created by Tom Lenaerts on 11/02/2022.
//

#ifndef ctpdata_hpp
#define ctpdata_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <map>

using namespace std;

class CtpEntry{
public:
    CtpEntry(){}
    CtpEntry(unsigned steps):_steps(steps){
        _probs.clear();
        for(unsigned i=0;i<=steps;i++){
            _probs.push_back(0.0);
        }
    }
    CtpEntry(unsigned steps, vector<double> data):_steps(steps){
        _probs.clear();
        for(unsigned i=0;i<data.size();i++){
            _probs.push_back(data[i]);
        }
    }

    CtpEntry(const CtpEntry& other);
    CtpEntry& operator=(const CtpEntry& other);

    ~CtpEntry(){
        _probs.clear();
    }

    unsigned steps() const {return _steps;}
    unsigned int size() const {return (unsigned)_probs.size();}
    double operator[](unsigned pos) const;
    void set(unsigned pos, double val);
    
    bool operator==(const CtpEntry& other) const;
    bool operator!=(const CtpEntry& other) const;

    std::vector<double>::iterator begin() {return _probs.begin();}
    std::vector<double>::iterator end() {return _probs.end();}
    std::vector<double>::const_iterator begin() const{return _probs.begin();}
    std::vector<double>::const_iterator end() const {return _probs.end();}

    friend std::ostream & operator<<(std::ostream &o, CtpEntry& s){return s.display(o);}
        
protected:
    virtual ostream& display(ostream& os) const ;

    unsigned _steps;
    vector<double> _probs;
};

class CtpData {
public:
    CtpData(){
        _dict.clear();
    }
    ~CtpData(){
        _dict.clear();
    }
    
    void add(string key, CtpEntry* value);
    CtpEntry* operator[](string key) const;

    friend std::ostream & operator<<(std::ostream &o, CtpData& s){return s.display(o);}

protected:
    virtual ostream& display(ostream& os) const ;
    map<string,CtpEntry*> _dict;
    
};
#endif /* ctpdata_hpp */
