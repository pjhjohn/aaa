#ifndef _LK_H
#define _LK_H

#include "macro.h"
#include "tsplib_io.h"
#include "tour.h"

	/// 1이면 속도가 빨라지는 대신에 질이 조금 떨어진다.
#define LK_MORE_SPEED        1
#define LK_DEPTH             40

	/////////////////////////////////////////////////////////////////////
	/** Lookbit Queue 구현한 클래스.
		Johnson paper에 나오는 lookbit queue를 구현하였다. 이 클래스는
		LK 알고리즘, 2-Opt, 3-Opt 등과 같은 TSP의 지역 최적화 알고리즘의
		수행 속도를 빠르게 하는 데 이용된다.
	*/
class CLookbitQueue
{
private:
		// lookbit 큐
	int*            m_look_bit_queue;
		// lookbit 큐의 크기
	int             m_size;
		// lookbit 큐의 머리 부분
	int             m_look_head;
		// lookbit 큐의 꼬리 부분
	int             m_look_tail;
	
public:
		// Lookbit Queue를 만든다. size는 도시 수.
	CLookbitQueue(int size);
	~CLookbitQueue();

		/* lookbit 큐에 초기 요소를 삽입한다.
		   tour의 도시 순서대로 요소를 삽입한다. tour안의 어떤 도시의 
		   incident edges가 c1,c2에 공통적으로 들어 있을 때는, 이 도시를
		   알고리즘에서 조사할 필요가 없다. 그러므로 이 도시는 큐에 
		   들어가지 않는다.
		   @param c1,c2 둘 다 지역 최적점 투어
		*/
	void    construct(C2EdgeTour* tour,
					  C2EdgeTour* c1=NULL, C2EdgeTour* c2=NULL);
		// 도시 c를 look bit 큐에 추가한다.
	void    addLookbit(int c);
		// 큐에서 도시 하나를 꺼낸다. 반환값은 도시 번호 또는 <0인 값
	int     deleteLookbit();
};

	////////////////////////////////////////////////////////////////////////
	/** 세그먼트 트리 클래스.
		Johnson의 paper에서 나온대로 구현하였다. order형식의 투어나
		C2EdgeTour같은 형식의 투어는 각각 flip과 between명령어의 시간 복잡도가
		$K\cdot n$이다($K$의 값은 크다). 이 때 세그먼트 트리를 사용하면
		K의 값을 몇십배로 줄일 수 있다.
	*/
class CSegmentTree
{
private:
		// segment tree가 다루고 있는 투어
	C2EdgeTour*     _tour; 
		// 투어의 도시 개수;
	int             _n;
	
		// 세그먼트를 나타내는 구조체
	struct SEGMENT
	{
			// reversal bit. 
		int         reversal;
			// 이 세그먼트의 첫번째 도시
		int         first;
			// 이 세그먼트의 마지막 도시
		int         last;
			// 트리로 간주될 때 왼쪽 자식
		SEGMENT*    left_child;
			// 트리로 간주될 때 오른쪽 자식
		SEGMENT*    right_child;
			// 리스트로 간주될 때 자기 앞의 세그먼트
		SEGMENT*    prev;
			// 리스트로 간주될 때 자기 뒤의 세그먼트
		SEGMENT*    next;
	};
		// 세그먼트 트리가 저장된다.
	SEGMENT         _seg_tree[(LK_DEPTH+2)*2];
		// 세그먼트 트리의 노드 수
	int             _seg_size;
		// 도시 번호 -> 순서가 나오는 배열.
	int*            _city_order;
		// 도시 번호 -> 이 도시가 포함된 세그먼트의 포인터
	SEGMENT**       _loc_pointer;
		// 도시 번호 -> #_loc_pointer[도시]#가 마지막으로 세팅된 시간. 
	int*            _time_stamp;
		// 현재 시간. #_loc_pointer[도시]#는 그 것이 마지막으로 세팅된 시간이 현재 시간과 같아야만 유효하다.
	int             _current_time;
	
public:
	CSegmentTree(int size);
	~CSegmentTree();

		// 시간 도장을 초기화한다.
	void        initTimeStamping();
		// 도시의 순서 번호를 구성한다.
	void        setupCityOrder(C2EdgeTour* tour);
		// t배열을 보고 i 레벨만큼 세그먼트 트리를 키운다. LK에서만 사용됨.
	void        setupTree(int i, int t[]);
		// 세그먼트 트리안에서 2-change
	void        do2Change(int a, int b, int c, int d);
		// 세그먼트 트리안에서 3-change
	void        do3Change(int t1, int t2, int t3, int t4,
							 int t5, int t6);
		// 세그먼트 트리안에서 4-change
	void        do4Change(int t1, int t2, int t3, int t4,
							 int t5, int t6, int t7, int t8);
		// city의 다음 순서에 있는 도시를 반환한다.
	int         getNext(int city);
		// city의 앞 순서에 있는 도시를 반환한다.
	int         getPrev(int city);
		// start, end도시 사이에 middle이 있는가?
	int         isBetween(int start, int middle, int end);
		// segment tree의 방향을 반대로 바꾼다.
	void        reverseTree();

private:
		// 도시 번호로 세그먼트를 찾는다.
	CSegmentTree::SEGMENT*findSegment(int city);
		// 도시의 순서 번호로 도시 번호를 찾는다.
	int         getCityByOrder(int order, int referCity);
		// 도시 c1, c2사이를 두 개의 세그먼트로 나눈다.
	void        splitSegment(int c1, int c2);
		// start, end사이의 세그먼트 리스트들을 뒤바꾼다.
	void        reverseSegment(SEGMENT* start, SEGMENT* end);
		// 제대로 된 트리인가? debug함수.
	int         isCorrect(); 
};

	/////////////////////////////////////////////////////////////////////////
	/** TSP를 풀기 위한 국지 최적화 알고리즘 구현을 위한 부모 클래스.
		LK(Lin-Kernighan Algorithm), 2-Opt, 3-Opt와 같은 알고리즘의
		공통 분모를 가지고 있다.
	*/
class CTSPLocalOpt
{
public:
		// 생성자. #num_city#는 도시의 수, #num_nn#은 조사할 이웃 도시의 수
	CTSPLocalOpt(int num_city, int num_nn);
	virtual ~CTSPLocalOpt() { }
	
		/* 국지 최적화 알고리즘을 수행하는 인터페이스 함수.
		   op1, op2, op3는 실제로 이 함수를 구현하는 국지 최적화 알고리즘마다
		   마음대로 사용할 수 있다. 즉 각각의 알고리즘마다 다른 의미로
		   사용할 수 있다.
		   @param tour                국지 최적화 알고리즘이 적용될 투어.
		   @param op1,op2,op3,op4,op5 옵션 투어들.
		*/
	virtual void    run(C2EdgeTour* tour, C2EdgeTour* op1=NULL,
						C2EdgeTour* op2=NULL, C2EdgeTour* op3=NULL,
						C2EdgeTour* op4=NULL, C2EdgeTour* op5=NULL) = 0;

protected:
		// Segment tree 변수
	CSegmentTree    m_segTree;
		// Lookbit 큐 변수
	CLookbitQueue   m_lookbitQueue;
		// 도시의 수
	int             m_numCity;
		// 이웃 도시의 수
	int             m_numNN;
};

	////////////////////////////////////////////////////////////////////////
	/// LK 알고리즘을 구현한 클래스
class CLK : public CTSPLocalOpt
{
public:
		/// 생성자. #num_city#는 도시 수, #num_nn#은 조사할 이웃 도시의 수
	CLK(int num_city, int num_nn); 
	~CLK();

		/* LK알고리즘을 수행하는 인터페이스 함수.
		   p1, p2, best는 전부 다 NULL이거나 전부 다 NULL이 아니어야 한다.
		   @param tour      LK알고리즘이 적용될 투어
		   @param p1,p2     tour를 만드는 데 이용한 부모들.     
		   @param best      지금까지 가장 좋은 해
		*/
	void    run(C2EdgeTour* tour, C2EdgeTour* p1=NULL,
				C2EdgeTour* p2=NULL, C2EdgeTour* best=NULL,
				C2EdgeTour* worst=NULL, C2EdgeTour* other=NULL);
	
private:
		/**@name LK알고리즘의 핵심 부분.
		   Lin과 Kernighan의 논문을 구현하였다.
		*/  //@{ start block
		/// runLK()에서 나온 tour
	C2EdgeTour*     m_tour;
		// 현재 레벨
	int             _i;
		// $G^*$를 얻기 위한 레벨
	int             _k;    
		// 현재 레벨까지의 이득(gain).
	double            *_G;     
		// $G^*$. 지금까지 가장 좋은 gain.
	double _Gstar; 
		// t 배열. $t_1, t_2, t_3, ... $
	int*            t;
	
		// 도시의 수
	int             _n;
		// 가장 가까운 이웃 도시의 수
	int             _nnn;  
	
        // $t_3$가 될 수 있는 후보를 저장하는 구조체
	struct LK_T3_CAND
	{
        int t3; int t4; int alter_t4; double gain;
	};
		// $t_3$ 후보의 개수
	int             _num_t3_cand_list;
		// $t_3$ 후보들이 저장되어 있는 배열
	LK_T3_CAND*     _t3_cand_list;
	
        /** $t_5$가 될 수 있는 후보를 저장하는 구조체.
			#get_t5_cand_list()#에서는 t5, t6, gain변수만 이용하고,
			#get_alter_t5_cand_list()#에서는 다 쓰인다. */
	struct LK_T5_CAND
	{
        int t5; int t6; int alter_t6; int t7; int t8;
        int code; double gain;
	};
		// $t_5$ 후보의 개수
	int             _num_t5_cand_list;
		// $t_5$ 후보들이 저장되어 있는 배열
	LK_T5_CAND*     _t5_cand_list;

		// 주어진 $t_1$으로 실제로 LK를 실행시킨다. improved는 향상 횟수
	double do_lk_search(int t1, int improved);
		// 더 깊은 레벨에서 수행한다. improved는 연속 향상 횟수
	void    search_deeper(int improved);
		// $G^*$를 바뀔 수 있는지 새로 평가한다.
	int     check_and_update(int improved);
		// 현재 레벨에서 k레벨로 투어를 되돌려 $G^*$의 이득을 얻도록 한다.
	void    reverse_change_to_best(int improved);
		// 현재 레벨에서 to레벨로 투어를 되돌린다.
	void    reverse_change(int to);
		// $t_3$후보들을 구한다.
	void    get_t3_cand_list();
		// $t_5$후보들을 구한다.
	void    get_t5_cand_list();
		// 대안 $t_5$후보들을 구한다.
	void    get_alter_t5_cand_list();

		// 2-change를 수행한다.
	void    make_two_change(int improved);
		// 3-change를 수행한다.
	void    make_three_change(int improved);
		// 4-change를 수행한다.
	void    make_four_change(int improved);
		//@} end block
};
#endif
	
