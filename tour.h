#ifndef _TOUR_H
#define _TOUR_H

	/** TSP(STSP, ATSP 둘 다)에서 사용되는 투어들의 부모 클래스. 
		STSP, ATSP 둘 다에서 공통적으로 필요한 함수들을 가지고 있다.
		ATSP의 경우 (v1, v2)와 (v2, v1)에지가 다르게 취급되므로,
		함수에서 두 도시를 인수로 받는 경우 순서를 생각하여야 한다. */
class CTour
{
 private:
		// findEdgeFirst, findEdgeNext에서 쓰이는 traversal을 위한 내부 변수.
	int        m_edge_start;
	int        m_edge_v1; 

protected:
		// 투어 비용
	double        m_length;
		// 도시의 수 
	int        m_size;
public:
		// 생성자
	CTour();

		/**@name 가상 함수들.
		   자식 클래스들은 이 함수들을 구현해야 한다. 특히, 순수 가상함수들은
		   반드시 구현되어야 한다. */ 
		//@{ start DOC++ block
		/// 이 투어와 tour사이의 해밍 거리를 반환한다.
	virtual int  getHammingDist(CTour* tour);
		// 랜덤 투어를 만든다.
	virtual void makeRandomTour()=0;
		// Edge(v1, v2)가 있는가?
	virtual int  isThereEdge(int v1, int v2)=0;
		// 도시를 나열하기 위하여 준비하는 함수. 도시 start부터 나열된다.
	virtual void findFirst(int start)=0;
		// 다음 도시를 반환한다.
	virtual int  findNext()=0;
		// 완전하고 적법한 투어인가?
	virtual int  isTour()=0;
		/** order형식으로 된 투어를 이 타입의 투어로 변환한다.
			@param src order형식의 투어
			@param size src의 크기. #m_size#와 같아야 한다.
		*/
	virtual int  convertFromOrder(int* src, int size)=0;

		// two-change를 하는 함수
	virtual void make2Change(int v1, int n1, int v2, int n2) {}
		// three-change를 하는 함수
	virtual void make3Change(int v1, int n1, int v2, int n2,
							 int v3, int n3) {}
		// four-change를 하는 함수
	virtual void make4Change(int v1, int n1, int v2, int n2,
							 int v3, int n3, int v4, int n4) {}
		//@} end of DOC++ block

		/// 투어 비용을 계산한다.
	double evaluate();
		// 도시의 개수를 리턴한다.
	int          getSize() { return m_size; };
		// 투어 비용을 반환한다.
	int          getLength() { return m_length; };
		// 에지를 나열하기 위해서 준비한다. 도시 start부터 시작한다.
	void         enumEdgeFirst(int start);
		/** 다음 에지를 반환한다.
			@param v1,v2 반환될 에지의 양쪽 도시가 들어간다. ATSP의 경우,
			에지는 v1->v2인 에지를 말한다.
			@return 성공할 때는 TRUE, 실패할 때는 FALSE
		*/
    int          enumEdgeNext(int* v1, int* v2);
		// 현재의 투어를 dest에 order형식으로 변환한다. size는 dest의 크기.
	void         convertToOrder(int* dest, int size);
		// 랜덤 투어를 dest에 order형식으로 변환한다. size는 dest의 크기.
	void         makeRandomOrderTour(int* dest, int size);
};

	/** STSP를 위한 투어 클래스.
		STSP 투어에서는 선행자나 계승자의 의미가 없다. 단지 모든 도시들이 2의
		degree를 가지면서 subcycle이 없으면 완전하고, 적법한 투어가 된다.
		그래서 이 클래스는 어떤 도시의 최대 연결된 도시의 수(degree)가 2임을
		이용하여, #m_e1, m_e2#배열로 투어를 표현한다. 

		#m_e1, m_e2#는 어떤 도시의 이웃도시들을 저장한다. 예를 들면,
		10 - 7 - 5 - ... 의 경우 #m_e1[7] = 10, m_e2[7] = 5, #
		#m_e2[10] = 7, m_e2[5] = 7, ...#처럼 표현 될 수 있다. 

		이 클래스의 초기화 방법은 2가지로 나눌 수 있다. 첫번째는 1단계 생성자를
		불러 한번에 초기화하는 방법이고 두번째는 2단계 생성자를 부르고 나서,
		create()함수를 부르는 방법이다.
	*/
class C2EdgeTour : public CTour
{
public:
		// 도시 i -> 도시 j. 도시 i와 j는 연결되어 있다. j<0이면 도시 i는 #m_e1#에서 다른 도시와 연결되어 있지 않다.
	int*        m_e1;
		// 도시 i -> 도시 k. 위와 마찬가지.
	int*        m_e2;
	
private:
		// findFirst, findNext에서 사용되는 내부 변수
	int         m_end;
	int         m_prev;
	int         m_cur; 
		// 멀티 그래프를 허용하는지를 나타내는 플래그. 
	int         m_isMult;

public:
		// 1단계 생성자
	C2EdgeTour(int size);
		// 2단계 생성자
	C2EdgeTour();
		// 투어의 구조를 생성한다. 
	void       create(int size);
	~C2EdgeTour();
	
		/**@name 가상 함수들
		   CTour의 순수 가상 함수들을 계승한 함수들이다. */
		//@{    block start
		///
	void        makeRandomTour();
	int         isThereEdge(int v1, int v2);
	void        findFirst(int start);
	int         findNext();
	int         isTour();
	int         convertFromOrder(int* src, int size);
	void        make2Change(int v1, int n1, int v2, int n2);
	void        make3Change(int v1, int n1, int v2, int n2,
							int v3, int n3);
    void        make4Change(int v1, int n1, int v2, int n2,
							int v3, int n3, int v4, int n4);
		//@}    block end

		/**@name C2EdgeTour의 순수한 연산 함수 */
		//@{    block start
		/// (v1, v2)에지를 투어에 더한다. 사전에 이 에지가 투어에 없어야 한다.
	void        addEdge(int v1, int v2);
		// (v1, v2)에지를 안전하게 삽입한다. 이미 에지가 있으면 No operation.
	void        addEdgeSafely(int v1, int v2);
		// (v1, v2)에지를 지운다.
	void        deleteEdge(int v1, int v2);
		// 모든 에지를 지워서, 투어가 비게 한다.
	void        deleteAllEdge();
		// 도시 v1의 degree를 리턴한다.
	int         getEdgeSize(int v1);
		// 투어에서 사용된 에지의 수를 리턴한다. 리턴값이 도시의 수와 같아야만 완전한 투어의 충분 조건이 된다.
	int         getNumEdge();
		// (v1, v2) 에지가 있는지를 빠르게 찾는다.
	int         isThereEdge2(int v1, int v2);
		/* 도시 v에 접해 있는 이웃 도시를 리턴한다.
		   @param v1,v2 리턴될 이웃 도시. 없으면 음수값이다. */
	void        getNeighbor(int v, int* v1, int* v2);
		// 도시 cur의 이전 도시를 prev라고 간주할 때, cur의 다음 도시를 리턴.
	int         getNext(int prev, int cur);
		// (v1, v2)에지가 (v3, v4)에지랑 같은 에지인가?
	int         is_same_edge(int v1, int v2, int v3, int v4);
		// 이 투어에 대한 연산도중 멀티그래프 형식을 허용할지 결정한다.
	void        setMultiGraphFlag(int flag);
		// src 투어를 이 투어에 복사한다.
	const C2EdgeTour& operator=(const C2EdgeTour& src);

		/* 도시 cur와 end사이에 city가 있는가?.
		   투어를 prev에서 cur로 가는 방향성 투어라고 가정할 때,
		   cur에서 출발하여 end까지 가는 패쓰안에 city가 있는 지의 여부를 
		   반환한다.
		   @param prev cur와 함께 투어의 방향을 정하는 데 쓰이는 도시
		   @param cur  prev와 함께 투어의 방향을 정한다. 그리고, 조사할 패쓰의
		   시작점이다.
		   @param city 패쓰안에 있는지를 알아 볼 도시 
		   @param end  패쓰의 끝점 */
	int         isBetween(int prev, int cur, int city, int end);
		// subcycle이 존재하는가?
	int         isDisjointCycle();
		//@}    block end

		/**@name 이 투어가 불완전한 투어일 때 사용되는 함수 */
		//@{    block start
		/** 불완전한 지금의 투어를 greedy방법으로 완전한 투어를 만든다.
			@param tabu tour의 금지 에지
			@param t1,t2 tour의 또다른 금지 에지. 옵션이다.
			@param useNN greedy를 빠르게 하기 위해서 이웃 도시만을 써라. option
		*/
		// 불완전한 투어를 임의의 연결방법으로 완전한 투어를 만든다.
	void        connectRandom(C2EdgeTour* tabu);
		// Nearest Neighbor방식으로 투어를 만든다. t1, t2는 또 다른 금지 옵션.
		/* 금지 에지를 만든다.
		   금지 에지란, 삽입했을 때 투어에 subcycle이 생기게 하는 에지를
		   말한다. 
		   @param tabu 만들어질 금지 에지 */
	void        constructTabuEdge(C2EdgeTour* tabu);
		// 에지(v1, v2)를 금지에지에 걸리지 않고, 이 투어에 삽입할 수 있는가?
	int         canAddEdge(C2EdgeTour* tabu, int v1, int v2);
		// 타겟 투어에 (v1, v2)을 삽입한다. 이 때 금지 에지의 내용도 바꾼다.
	void        addEdgeWithTabu(C2EdgeTour* tabu, int v1, int v2);
		/* 도시 v의 끝점을 리턴한다.
		   투어가 불완전할 때, 투어는 여러 개의, 구별되는 패스로 이루어져 있다.
		   이 때, 도시 v가 한 패스의 끝점이라면 이 함수는 그 패스의 다른
		   끝점을 리턴한다. v가 패스의 중간이라면 음수값을 리턴한다.
		   v의 degree가 0일때도 음수값을 리턴한다. */
	int         getEndPoint(int v);
		//@}    block end
	
};
	
inline void C2EdgeTour::deleteAllEdge()
{
	for( int i=0; i<m_size; i++) m_e1[i] = m_e2[i] = -1;
}

inline void C2EdgeTour::getNeighbor(int v, int* v1, int* v2)
{
	*v1 = m_e1[v];	    *v2 = m_e2[v];
}

inline int  C2EdgeTour::getNext(int prev, int cur)
{
	return (m_e1[cur] == prev) ? m_e2[cur] : m_e1[cur];
}


inline int  C2EdgeTour::is_same_edge(int v1, int v2, int v3, int v4)
{
	return ((v1==v3&&v2==v4) || (v1==v4&&v2==v3));
}

inline int  C2EdgeTour::isThereEdge2(int v1, int v2)
{
    return ((m_e1[v1] == v2) || (m_e2[v1] == v2));
}

inline void C2EdgeTour::setMultiGraphFlag(int flag)
{
	m_isMult = flag;
}

inline int  C2EdgeTour::getEdgeSize(int v1)
{
	if( m_e1[v1]>=0 && m_e2[v1]>=0)      return 2;
	else if( m_e1[v1]>=0 || m_e2[v1]>=0) return 1;
	else                                 return 0;
}

#endif













