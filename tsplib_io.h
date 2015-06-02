#ifndef TSPLIB_IO_H
#define TSPLIB_IO_H

	/**@name TSPLIB 입출력 루틴.
	   \URL[TSPLIB]{ftp://softlib.rice.edu/pub/tsplib} 파일을 읽어 TSP
	   인스턴스를 구성하고, 그 것을 접근하는 데 이용되는 변수와 함수 집합을
	   제공한다.
	   
	   밑에서 기술된 order형식의 tour란, 예를 들면 7->5->3->... 이라는 투어를
	   tour[0] = 7, tour[1] = 5, tour[2] = 3, ... 식으로 구성된 투어를 말한다.
	*/
	//@{  Start of DOC++ block

#ifndef DISTANCE_CALC_METHOD
	/** 거리 데이터를 저장하는 방법을 지정하는 컴파일러 플래그. 
		이 매크로의 값에 따라 거리 데이터를 접근하는 방법이 달라진다. 도시의
		개수가 커질수록 거리를 저장하는 배열도 커지게 되어 10000개이상의
		도시에서는 더이상 Matrix로 저장할 수 없기 때문에 다른 방법을 생각해야
		한다. 이 매크로는 거리를 저장하거나 계산하는 방법을 지정하는 것이다.
		이 매크로는 내부적으로만 사용되어야 한다. 즉 #tsplib_io.cc#에서만
		사용되는 것이다.
		\begin{itemize}
		\item 0이면, 거리를 Full Matrix에 저장한다. 즉 이 행렬은
		$4n^2$byte를 차지한다. 계산이 간단하므로 접근 속도는 제일 빠르다.
		\item 1이면, 거리를 Half Matrix에 저장한다. 이 방법은 symmetric TSP에만
		적용된다. $\frac{4(n-1)(n-2)}{2}$byte를 차지한다. 즉 Full Matrix보다
		공간을 절반정도 절약한다. 속도도 빠른 편이다.
		\item 2이면, 거리를 온라인으로 계산한다. 즉 거리가 필요할 때마다
		주어진 도시의 좌표로 실수 계산을 하기 때문에 정말 느리다. 도시의 개수가
		클 때에만 이 방법을 써야 한다. 실제로는 거리를 저장하는 데 캐시를
		사용하기 때문에 약 $4((neighbors+\alpha)*n)$byte만큼 차지한다. 이
		방법으로 계산할 때는 dist()매크로는 ConstructNN()함수 호출후부터
		유효하다. 나머지 방법은 ReadTspFile()함수 호출후부터 유효하다.
		\end{itemize}
	*/
#define DISTANCE_CALC_METHOD    0
#endif


#if     (DISTANCE_CALC_METHOD==0)
	/// 도시 c1, c2사이의 거리를 반환한다.
#define dist(c1, c2)        gDistMat[(c1)*gNumCity + (c2)]
#elif   (DISTANCE_CALC_METHOD==1)
#define dist(c1, c2)		gDistMat[max(c1,c2)*(max(c1,c2)-1)/2 + min(c1,c2)]
#elif   (DISTANCE_CALC_METHOD==2)
#define dist(c1, c2)        calcOnline(c1, c2)
#endif

	/// 도시 i에서 nth번째 떨어진 도시의 번호
#define nni(i, nth)         gNNI[(i)*gNumNN + (nth)]
	/// 도시 i에서 nth번째 떨어진 도시와의 거리
#define nnd(i, nth)         (dist((i), nni(i, nth)))

	/// 문제의 인스턴스가 symmetric TSP임을 가리킨다.
#define TYPE_STSP           0
	/// 문제의 인스턴스가 asymmetric TSP임을 가리킨다.
#define TYPE_ATSP           1

	// 도시의 수
extern unsigned int    gNumCity;
	// 인스턴스의 타입 : symmetric 또는 asymmetric TSP인지 지정
extern int             gType; 
	// 저장된 가장 가까운 이웃 도시의 수
extern int             gNumNN; 
	// 거리 행렬. 직접 사용하지 말고, dist매크로를 이용한다.
extern double*            gDistMat; 
	// 가장 가까운 이웃 도시 번호들을 저장한 배열. nni매크로를 이용하라.
extern int*            gNNI;           

	/** TSP파일을 읽고 인스턴스를 구성한다.
		@param graphname 인스턴스 파일이름. 확장자는 .tsp 또는 .atsp가 되어야
		한다.  */

void ReadTspFile(char *graphname);
	/** 가장 가까운 이웃 도시들의 배열을 만든다.
		@param numNN 한 도시에 대해서 만들어질 이웃 도시들의 수
		@param is_quad 이웃 도시들을 계산할 때 도시들을 4개의 구역으로 나누어서
		계산할 지 여부를 결정한다. 인스턴스가 2차원일때만 가능하다.
	*/
void ConstructNN(int numNN, int is_quad=0);

	// 인스턴스 파일에 대한 정보를 출력한다.
void PrintTspInfo();

	/** 인스턴스 파일이 n차원일 때, 투어의 순열대신에 좌표의 순열을 출력한다.
		@param fd 파일 기술자
		@param tour order형식으로 된 여행 경로
		@param size tour의 길이	 */
void PrintCoords(FILE* fd, int* tour, int size);

	/** 인스턴스 파일이 n차원일 때, 한 도시의 좌표를 얻어 낸다.
		@param city 도시 번호
		@param pt   반환될 좌표가 들어갈 배열
		@param size pt배열의 크기. n차원보다 더 커야 한다. 
		@return 차원의 수를 리턴한다. 즉 n을 리턴하고 
		실패하면, -n을 리턴한다.  */
int GetCoords(int city, double* pt, int size);

	/** 투어 파일을 읽는다.
		@param file 투어 파일 이름
		@param tour 반환될 투어가 order형식으로 들어갈 배열. 
		@param size tour의 크기.
		@return 성공 여부를 알리는 boolean값 */
int ReadTourFile(char* file, int *tour, int size);

	/** 최적의 투어 파일을 읽는다.
		@param opt_tour 반환될 투어가 order형식으로 들어갈 배열. 
		@param size opt\_tour의 크기. 도시 수보다 크거나 같아야 한다.
		@return 성공 여부를 알리는 boolean값 */
int ReadOptTourFile(int *opt_tour, int size);

	// 인스턴스의 최적 비용을 반환한다. 모를 경우에는 0을 반환한다.
int GetOptimumCost();

	/** 도시 번호를 재구성한다. 
		@param remap 리맵된 도시 -> 원래의 도시 번호를 돌려주는 배열이다. */
void RemapCities(int* remap);

	// 리맵된 도시의 원래 도시 번호를 반환한다.
int GetOrgCity(int city);

	/** 투어를 출력한다. 만약 도시들이 리맵되어 있다면, 원래의 도시 번호들로
		변환하여 출력한다.
		@param fd 파일 기술자
		@param tour 출력할 order형식의 투어 배열.
		@param size tour의 크기
		@param name 투어 내용에 들어 갈 파일 이름
		@param comment 주석문
	*/
void PrintTour(FILE* fd, int* tour, int size, char* name=NULL,
			   char* comment=NULL);

	// 인스턴스의 그래프 이름을 dest에 복사한다.
void CopyGraphName(char* dest);

	// 도시 i와 j의 사이의 거리를 반환한다. dist매크로를 사용할 것. 
int calcOnline(int i, int j);

	//@}  End of DOC++ block
#endif /* end of file */









