/*
 * =====================================================================================
 *
 *       Filename:  ga_solver.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24/11/2011 12:27:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eduardo Feo
 *        Company:  IDSIA
 *
 * =====================================================================================
 */


#ifndef _GA_SOLVER_H
#define _GA_SOLVER_H

#include "main.h"
#include <ga/ga.h> // and the 2D binary string genome

#include "solver.h"

class Network;
class Graph;
class Matheuristic;
class MyCplexSolver;


class GASol
{
  double m_cpuTimeSolving;
  bool m_hasSolution;
  RNPSolutionPtr m_sol;
  Matheuristic *m_matheuristic;
  //vector< pair< double, double> > objectiveTimeline;
  public:
  GASol(Network *,Network *);
  double solve();
  bool hasSolution(){ return m_hasSolution;}
  double test();
  RNPSolutionPtr solution(){ return m_sol;}
  void testRandom();
  void setMatheuristic(Matheuristic *mth){ m_matheuristic = mth;}
  /*void addObjectiveTime(double time, double obj)
    { 
    objectiveTimeline.push_back(make_pair(time, obj));
    }
    vector< pair<double, double> > getObjectiveTimeline(){ return objectiveTimeline;}
    */

  //Graph* create_graph_from_sol();

  //void write_time_to_file(string filename);
  //void write_solution_to_file(string filename);
};


class MyGenome : public GA1DArrayGenome< int >
{
  public:

    GADefineIdentity("MyGenome", 201);

    static void Init(GAGenome&);
    static int Mutate(GAGenome&, float);
    static float Compare(const GAGenome&, const GAGenome&);
    static float Evaluate(GAGenome&);
    static int Cross(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
    static int repairs;

    string toString(string separator="");

    MyGenome(int k) : GA1DArrayGenome<int>(k+1) { 
      initializer(Init);
      mutator(Mutate);
      comparator(Compare);
      evaluator(Evaluate); 
      crossover(Cross); 
    }
    void fromLpSolution(LpSolutionPtr lpsol);


    void repair();

    MyGenome(const MyGenome& orig): GA1DArrayGenome<int>(orig.length())  { 
      MyGenome::copy(orig); 
      //printf("Copied and evaluated %d\n", _evaluated);
    }
    virtual ~MyGenome() {

    }
    virtual int write (STD_OSTREAM & os) const {
      for(unsigned int i=1; i<=gene(0); i++)
	os << gene(i) << " ";
      return 0;
    }
    MyGenome& operator=(const GAGenome& orig){
      if(&orig != this) MyGenome::copy(orig);
      return *this;
    }

    virtual GAGenome* clone(CloneMethod c) const {
      //printf("Cloning and evaluated %d\n",_evaluated);
      MyGenome *cloned = new MyGenome(*this);
      return cloned;
    }
    virtual int equal(const GAGenome&) const;
    virtual void copy(const GAGenome& orig) {
      MyGenome &morig = (MyGenome&)orig;
      GA1DArrayGenome<int>::copy(morig);  // this copies all of the base genome parts
      // copy any parts of MyObject here
      // copy any parts of MyGenome here
    }

    // any data/member functions specific to this new class
};


#endif
