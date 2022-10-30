//
//  ctpdata.cpp
//  ctpbeliefs
//
//  Created by Tom Lenaerts on 11/02/2022.
//

#include "ctpdata.hpp"
#include <cmath>

CtpEntry::CtpEntry(const CtpEntry& other){
    _steps=other.steps();
    _probs.clear();
    for(unsigned i=0; i < other.size(); i++){
        _probs.push_back(other[i]);
    }
}

CtpEntry& CtpEntry::operator=(const CtpEntry& other){
    _steps=other.steps();
    _probs.clear();
    for(unsigned i=0; i < other.size(); i++){
        _probs.push_back(other[i]);
    }
    return *this;
}


double CtpEntry::operator[](unsigned pos) const{
    if(pos>=0 && pos < _probs.size())
        return _probs[pos];
    return NAN;
}

void CtpEntry::set(unsigned pos, double val){
    if(pos>=0 && pos < _probs.size()){
        _probs[pos] = val;
    }
}

ostream& CtpEntry::display(ostream& os) const {
    os << "(" << _steps << "){";
    for(unsigned i=0; i < _probs.size(); i++){
        os << _probs[i];
        if(i < (_probs.size()-1))
            os <<",";
    }
    os << "}";
    return os;
}

bool CtpEntry::operator==(const CtpEntry& other) const{
    bool ret=(_steps == other.steps() && _probs.size() == other.size());
    if(!ret)
        return ret;
    
    for(unsigned i=0; i < _probs.size(); i++){
        ret = ret && (_probs[i] == other[i]);
        if(!ret)
            return ret;
    }
    return ret;
}

bool CtpEntry::operator!=(const CtpEntry& other) const{
    return !(*this == other);
}


void CtpData::add(string key, CtpEntry* value){
    _dict[key]=value;
}

CtpEntry* CtpData::operator[](string key) const{
    map<string,CtpEntry*>::const_iterator elm = _dict.find(key);
    CtpEntry* ret=NULL;
    if(elm != _dict.end())
        ret = elm->second;
    return ret;
}


ostream& CtpData::display(ostream& os) const {
    os << "DICT=["<< endl;
    map<string,CtpEntry*>::const_iterator start = _dict.begin();
    map<string,CtpEntry*>::const_iterator stop = _dict.end();
    while (start != stop) {
        string tmp = start->first;
        CtpEntry* e = start->second;
        os << "\t["<<tmp<<","<<*e<<"]";
        os << endl;
        start++;
    }
    os << "]" << endl;
    return os;
}



