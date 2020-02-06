#ifndef PI
#define PI acos(-1)
#endif


#ifndef RANDOMNUMBER_H
#define RANDOMNUMBER_H

#include<cmath>
#include<random>
#include<map>
using std::map;

class RandomNumber{
public:
    RandomNumber(){
        rng.seed(std::random_device()());
    }
    virtual ~RandomNumber(){ }
    int getRandomInt(int down, int up){
        std::uniform_int_distribution<int> iequaldis(down, up);
        return iequaldis(rng);
    }
    double getRandomfloat(double down, double up){
        std::uniform_real_distribution<double> fequaldis(down, up);
        return fequaldis(rng);
    }
    double getRandomZtoO(){
        std::uniform_real_distribution<double> fequaldis(0, 1);
        return fequaldis(rng);
    }
    double getRandomfromDistribution(map<double,double> const& rhomap){
        std::uniform_real_distribution<double> fequaldis(0, 1);
        return rhomap.lower_bound(fequaldis(rng))->second;
    }

private:
    std::mt19937 rng;
};



#endif //RANDOMNUMBER_H