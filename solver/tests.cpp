#include "CpuTime.h"
#include "main.h"




int main(int argc, char *argv[])
{

	CpuTime cputime;

	cputime.start();

	
	for(int j=1000; j< 5000; j++)
	{
//		printf("%d\n",j);
		double x=j;

		std::vector<double> myvec;
		for(int i=0; i< 100000; i++)
		{
			double r = rand()*x;
			myvec.push_back(r);

		}
		std::sort(myvec.begin(), myvec.end());
	}
	double telapsed = cputime.cpu_time_elapsed();
	printf("Time elapsed %.2f\n", telapsed);

	return 0;
}



