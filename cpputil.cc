    /* Test driver for cpputil test :
       gcc cpputil.cc -DCPP_TEST -lm
       usage a.out tsp-file
    */

#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include "cpputil.h"

#define affirm(x)
#define trace

	//#define affirm(x)    assert(x)
	//#define trace        printf

	/*******************************************************************/
	/************************ CTimer class *****************************/
void CTimer::start()
{
    m_elapsedStart = times(&m_start);
}

float CTimer::stop(int want)
{
	float       interval, clk_tick;

	m_elapsedEnd = times(&m_end);
	clk_tick     = sysconf(_SC_CLK_TCK);
	if( want == 0) 
		interval = (float)(m_end.tms_utime-m_start.tms_utime)/clk_tick;
	else
		interval = (float)(m_elapsedEnd - m_elapsedStart)/clk_tick;

	return interval;
}

float CTimer::report(int want)
{
	float       interval, clk_tick;
	
	clk_tick     = sysconf(_SC_CLK_TCK);
	if( want == 0) 
		interval = (float)(m_end.tms_utime-m_start.tms_utime)/clk_tick;
	else
		interval = (float)(m_elapsedEnd - m_elapsedStart)/clk_tick;
	
	return interval;
}

	/****************************************************************/
	/************************ CStat class ***************************/
CStat::CStat(int mode)
{
	m_mode = mode;
	setup();
}

void CStat::setup()
{
	m_count = 0;
	m_sum = m_square_sum = 0;
	m_vavg = m_square_vavg = 0;
	m_min_value = INT_MAX*4.0;
	m_max_value = INT_MIN*4.0;
}

int CStat::add(double value)
{
	int res = 0;

	m_count++;
	if( value < m_min_value)
	{
		m_min_value = value; res += MIN_UPDATE; affirm( res & MIN_UPDATE);
	}
	if( value > m_max_value) 
	{
		m_max_value = value; res += MAX_UPDATE; affirm( res & MAX_UPDATE);
	}

	if( m_mode==0)
	{
		m_sum        += value;
		m_square_sum += (value*value);
	}
	else
	{
		affirm( m_count > 0);
		m_sum        += value;
		m_vavg        = (m_vavg*(m_count-1.0) + value)/((double)m_count);
		m_square_vavg = (m_square_vavg*(m_count-1.0) + value*value)/
			((double)m_count);
	}
	return res;
}

double CStat::getAvg()
{
	if( m_mode == 0)
		return (m_sum / (double)m_count);
	else
		return m_vavg;
}

double CStat::getStd()
{
	double avg = getAvg();

	if( m_mode == 0)
		return sqrt(fabs(m_square_sum/(double)m_count - avg*avg));
	else
		return sqrt(fabs(m_square_vavg - avg*avg));
}

double CStat::getCiv95()
{
	return (1.96*getGrpStd());
}

double CStat::getCiv99()
{
	return (2.58*getGrpStd());
}

double CStat::getGrpStd()
{
	double avg, std;

	avg = getAvg();	std = getStd();
	return std/sqrt(m_count);
}

double CStat::getAbovePercent(double base)
{
	double avg;
	if( base <= 0.0005)
		return -1;
	else
	{
		avg = getAvg();
		return (((avg - base)*100.0)/base);
	}
}

////////////////////// TEST DRIVER //////////////////////////////
#ifdef CPP_TEST
int main(int argc, char* argv[])
{
	CStat  stater;
	int    num;
	double datum;
	FILE*  fp;

	if( argc >= 2)
	{
		fp = fopen(argv[1], "r");
		
		while( (num = fscanf(fp, "%lf", &datum)) == 1)
		{
			stater.add(datum);
			printf(" %d : std = %lf\n", stater.getCount(), stater.getStd());
		}
		fclose(fp);
	}
	return 0;
}
#endif

