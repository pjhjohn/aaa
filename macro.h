#ifndef _MACRO_H_
#define _MACRO_H_

	/**@name 일반적으로 쓰이는 매크로들 */
	//@{ block start
extern "C"
{
#ifdef __sparc
    extern long int random(void);
	extern int gethostname(char *name, int namelen); 
    extern void srandom(unsigned int);
#endif
}


#ifndef SWAP
	/// x와 y의 값을 맞바꾼다.
#define SWAP(x, y, t)   {t=x; x=y; y=t;}
#endif 

	/// TRUE
#define TRUE            1
	/// FALSE
#define FALSE           0
	/// BOOL type
typedef int BOOL;

	/// 실수 x를 반올림하여 정수로 만든다.
#define nint(x)             ((int)((double)(x)+0.5))
	/// 실수 x의 소수점을 잘라서 정수로 만든다.
#define cint(x)             ((int)(ceil(x)))

	/// x와 y중에서 작은 값
#define min(x, y)           ((x)<(y)?(x):(y))
	/// x와 y중에서 큰 값
#define max(x, y)           ((x)<(y)?(y):(x))
		
	/// 정수쌍을 나타내는 범용 구조체
struct IntPair { int i1; int i2; };

	/// 정수 3중쌍을 나타내는 범용 구조체 
struct IntTriple { int i1; int i2; int i3; }; 

	/// (정수, 실수)를 나타내는 범용 구조체
struct IntFloat { int i; float f; };

	/// 실수쌍을 나타내는 범용 구조체
struct FloatPair { float f1; float f2; };

	/** 선택 정렬을 한다.
		적은 수의 데이터들을 소트할 때 이용하는 매크로 함수이다. increasing
		order 정렬과 decreasing order 정렬을 둘 다 할 수 있다. 다음은
		increasing order(정렬 후, 매처음 요소는 가장 작은 요소이다.) 정렬을
		하는 예이다.\\
		##\\
		#     ELEMENT data[20];#\\
		#     data[0].dist = 5; data[1].dist = 3; data[2].dist = 9; ...#\\
		#     SORT(ELEMENT, data, 20, data[cur].dist < data[key].dist); #\\
		##\\
		cur변수는 현재 비교되는 요소를 가리키는 인덱스이고, key변수는 가장
		크거나 작은 요소를 가리키는 인덱스이다. 
		위의 경우, key는 가장 작은 요소를 가리키고 있다.
		cur, key라는 변수를 이 매크로 주변에서 선언하면, side effect가 생길수
		있으므로 조심해야 한다.
		@param type         요소의 형
		@param arr_name     정렬될 요소가 있는 배열
		@param size         요소 개수
		@param replace_cond 현재 요소와 key요소를 맞바꾸게 되어야 하는 조건
	*/
#define SORT(type, arr_name, size, replace_cond)                \
		{                                                       \
			int i, cur, key;                                    \
			type temp_;                                         \
			for( i=0; i<size-1; i++)                            \
			{                                                   \
				key = i;                                        \
				for( cur=i+1; cur<size; cur++)                  \
					if( replace_cond) key = cur;                \
				SWAP(arr_name[i], arr_name[key], temp_);        \
			}                                                   \
		}

	/** 최대값 또는 최소값을 가진 배열의 색인을 알아낸다.
		어떤 구조체 또는 기본 타입의 배열중에서 조건에 따라 최대값 또는
		최소값을 가지는 요소의 색인을 지정된 변수에 복사한다.
		다음은 크기가 20인 data배열안에서 data[5~10]사이에 있는 요소 중
		최소값을 찾는 예이다(data[5], data[10]도 비교된다.) 
		minDistance라는 변수는 data배열중에서
		가장 작은 값을 가지는 요소의 색인값을 가지게 된다.
		##\\
		#     ELEMENT data[20];#\\
		#     data[0].dist = 5; data[1].dist = 3; data[2].dist = 9; ...#\\
		#     GET_MINMAX(minDistance, 5, 10, data[cur].dist<data[key].dist);#\\
		##\\
		cur변수는 현재 비교되는 요소를 가리키는 인덱스이고, key변수는 가장
		크거나 작은 요소를 가리키는 인덱스이다. 
		위의 경우, key는 가장 작은 요소를 가리키고 있다.
		cur, key라는 변수를 이 매크로 주변에서 선언하면, side effect가 생길수
		있으므로 조심해야 한다.
		@param dest         결과가 복사될 변수. int형이어야 한다.
		@param start,end    배열에서 조사될 범위의 시작과 끝 색인값
		@param cond         현재 요소와 key요소를 맞바꾸게 되어야 하는 조건
	*/
#define GET_MINMAX(dest, start, end, cond)                      \
        {                                                       \
			int cur, key;                                       \
			key = (start);                                      \
			for( cur=key+1; cur<=(end); cur++)                  \
				if( cond) key = cur;                            \
			dest = key;                                         \
		}

			//@} block end
#endif







