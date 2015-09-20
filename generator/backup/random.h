/*
 * =====================================================================================
 *
 *       Filename:  random.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/01/2010 01:03:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "generator.h"

template< typename TBoostDMode >
struct randgen {
	typedef TBoostDMode type;
	typedef std::vector< typename TBoostDMode::result_type > collection;

//	! Constructor
	randgen( typename TBoostDMode::input_type const & min,
		 typename TBoostDMode::input_type const & max )
		: gm_( min, max ),
		vg_(seed_, gm_) {	}

	// Retrieve one random number
	typename TBoostDMode::result_type operator()() { 
		return vg_();	
	}

	// Retrieves a vector of random numbers
	collection giveMe( int n )  {
		collection v( n );
		std::generate( v.begin(), v.end(), *this );
		return v;
	}

	private:
	static boost::mt19937 seed_;
	TBoostDMode gm_;
	boost::variate_generator< boost::mt19937&, TBoostDMode > vg_;
};

template<class TBoostDMode>
boost::mt19937 randgen<TBoostDMode>::seed_( std::time(0) );






