#ifndef _CPPUTIL_H_
#define _CPPUTIL_H_

#include <sys/types.h>
#include <sys/times.h>

	/** @name 작지만 유용한 클래스 */ 
	//@{            Start of DOC++ block
	/** 걸린 시간을 보고하는 클래스. */
class CTimer
{
private:
		// 기록을 시작한 시간
	struct tms  m_start;
		// 기록을 시작할 때에, times()의 리턴값
	clock_t     m_elapsedStart;
		// 마지막으로 기록을 중단한 시간
	struct tms  m_end;
		// 기록을 중단할 때에, times()의 리턴값
	clock_t     m_elapsedEnd;
public:
		// 시간 기록 시작
	void        start();
		// start()후 걸린 시간을 반환. which가 0이면 사용자 시간을, 1이면 실제 걸린 시간을 반환한다.
	float       stop(int which=0);
		// 마지막 stop()가 반환한 값을 그대로 돌려줌. which가 0이면 사용자 시간을, 1이면 실제 걸린 시간을 반환한다.
	float       report(int which=0);
};

#define MIN_UPDATE      1
#define MAX_UPDATE      2

	// 많은 수치 자료들의 평균따위를 내주는 클래스
class CStat
{
private:
		// 모드. 0이 아니면 가평균방식(더 큰 숫자 취급가능). 0이면 합산방식.
	int         m_mode;
		// 자료 수
	int         m_count;
		// 자료들의 합. $\sum x_i$. 합산 방식에서 사용됨. 
	double       m_sum;
		// 자료들의 제곱의 합. $\sum x_i ^2$. 합산 방식에서 사용됨.
	double       m_square_sum;
		// 자료들의 가평균(virtual avg). $\frac{\sum x_i}{n}$. 가평균 방식.
	double       m_vavg;
		// 자료들의 제곱의 가평균. $\frac{\sum {x_i}^2}{n}$. 가평균 방식.
	double       m_square_vavg;
		// 가장 작은 자료의 값
	double       m_min_value;
		// 가장 큰 자료의 값
	double       m_max_value;
public:
		// 생성자
	CStat(int mode=0);
		// 모든 변수들을 초기화한다
	void        setup();
		/* 데이터 하나를 추가한다. 
		   @param  value  추가할 데이터
		   @return <0-3>. 이 자료가 가장 큰 또는 작은 자료인지를 알려준다.
		   반환값과 MIN\_UPDATE, MAX\_UPDATE와 bitwise and를 하면
		   그 여부를 알 수 있다. */
	int         add(double value);
		// 자료수를 반환
	int         getCount() { return m_count; }
		// 가장 작은 자료 값
	double       getMin() { return m_min_value; }
		// 가장 큰 자료 값
	double       getMax() { return m_max_value; }
		// 자료들의 합을 반환
	double       getSum() { return m_sum; }
		// 자료들의 평균을 반환
	double       getAvg();
		// 표준 편차를 반환
	double       getStd();
		// 95\% 신뢰도를 가질때, 차이를 리턴
	double       getCiv95();
		// 99\% 신뢰도를 가질때, 차이를 리턴
	double       getCiv99();
		// 그룹 편차를 반환.
	double       getGrpStd();
		// 기준값에 비해 평균이 얼마나 높은지 백분율 값을 반환한다.
	double       getAbovePercent(double base);
};

	//@}       End of DOC++ block 
#endif























































