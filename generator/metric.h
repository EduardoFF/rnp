/*
 * =====================================================================================
 *
 *       Filename:  metric.h
 *
 *    Description:  Calculate metric for network topology
 *
 *        Version:  1.0
 *        Created:  11/12/2010 09:46:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eduardo Feo
 *        Company:  IDSIA
 *
 * =====================================================================================
 */


#ifndef _METRIC_H
#define _METRIC_H

class Metric
{
  public:
    Network *net;
    Metric(Network *n):net(n){}
    virtual double operator()(string i, string j)=0;
    virtual double average() const = 0 ;
    virtual double max() const = 0;
    virtual ~Metric(){}
    virtual void stats(ostream &os) const;
};

ostream &
operator<<(ostream &os, const Metric &m);

class Constant_Metric: public Metric
{
  public:
    double m_val;
    Constant_Metric(Network *n, double val):
      Metric(n),m_val(val){}
    double operator()(string i, string j){ return m_val;}
    double average() const { return m_val;}
    double max() const { return m_val;}
    ~Constant_Metric() {}
};


class EstPRR_Metric: public Metric
{
  public:
    double max_cost;
    map<pair<string, string>, double> metric;
    EstPRR_Metric(Network *n);
    double operator()(string i, string j);
    double average() const;
    double max() const;
    ~EstPRR_Metric() {}
};


class DistSQ_Metric: public Metric
{
  public:
    double max_cost;
    map<pair<string, string>, double> metric;
    DistSQ_Metric(Network *n);
    double operator()(string i, string j);
    double average() const;
    double max() const;
    ~DistSQ_Metric() {}
};


/**
 * \brief  uses weight loaded from network file
 */
class LinkWeight_Metric: public Metric
{
  public:
    double max_cost;
    map<pair<string, string>, double> metric;
    LinkWeight_Metric(Network *n);
    double operator()(string i, string j);
    double average() const;
    double max() const;
    ~LinkWeight_Metric() {}
};

#endif


