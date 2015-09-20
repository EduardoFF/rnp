
#include "main.h"

#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sys/stat.h>
#include <sys/wait.h>


#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>


#include "solver.h"
using namespace std;
void testBoostCompression()
{

	// open the archive
	std::ostringstream mystream;

	namespace bio=boost::iostreams;
	bio::filtering_stream<bio::output> f;

	f.push(bio::gzip_compressor());
	f.push(mystream);

	f << "Hello";
	f.flush();
}

void
testBoostSerialize()
{
  LpSolution *lpsol = new LpSolution();
  lpsol->value["TEST"] = 12345;
  lpsol->objval = 12345;
  int buf_len = 65536;
  int msg_size;
  std::vector<char> buff_vec;

  printf("Attempting to serialize lpsol %f\n", lpsol->objval);
  save_lpsolution(*lpsol, buff_vec, buf_len);
  printf("Serialized BINARY size %d\n", buff_vec.size());

}

int main()
{
	printf("Testing serialization\n");
	testBoostSerialize();
	return 1;
}
