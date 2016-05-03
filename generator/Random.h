/*
 * =====================================================================================
 *
 *       Filename:  Random.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/15/2010 04:40:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef _RANDOM
#define _RANDOM


class RandomGenerator
{
 public:
  virtual double operator()() = 0;
};

class UniformRandom : public RandomGenerator{
 public:
  variate_generator<mt19937, uniform_real<> > gen;

 UniformRandom(int seed) : gen(mt19937(seed),uniform_real<>(0,1)){};
  double operator()();
};

#endif
