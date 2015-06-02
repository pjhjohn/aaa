/* Test driver for LK test :
  gcc lk.cc edgetour.cc tsplib_io.cc cpputil.cc -DLK_TEST -lm -DDISTANCE_CALC_METHOD=1
  usage a.out tsp-file
*/
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

const double EPS = 1e-8;

#include "lk.h"

	//////////////////////////////// DEBUG ///////////////////////////////
#define affirm(x)
#define trace

//#define affirm(x)           assert(x)
//#define trace               printf

	//////////////////////// CTSPLocalOpt class ///////////////////////////
CTSPLocalOpt::CTSPLocalOpt(int num_city, int num_nn) 
	: m_lookbitQueue(num_city),	m_segTree(num_city)
{
	m_numCity = num_city;       m_numNN = num_nn;
}

	/*********************************************************************/
	/********************* LK class  *************************************/
	/********************* Interface Functions ***************************/
CLK::CLK(int num_city, int num_nn) : CTSPLocalOpt(num_city, num_nn)
{
    trace("Entering CLK()\n");
    _n                  = num_city; _nnn = num_nn;
    _G                  = new double[_n+1];
    t                   = new int[_n*2+2];
    _t3_cand_list       = new LK_T3_CAND[_nnn];
    _t5_cand_list       = new LK_T5_CAND[_nnn];

    trace("Quitting CLK()\n");
}

CLK::~CLK()
{
	trace("Entering ~CLK()\n");
    delete[] _G;
    delete[] t;
    delete[] _t3_cand_list;
    delete[] _t5_cand_list;

    trace("Quitting ~CLK()\n");
}

void CLK::run(C2EdgeTour* tour, C2EdgeTour* p1, C2EdgeTour* p2,
			  C2EdgeTour* best, C2EdgeTour* worst, C2EdgeTour* other)
{
	trace("Entering CLK::run\n");
    int     t1, improved;
	
    m_lookbitQueue.construct(tour, p1, p2);

	m_tour = tour;
	
		/* initialize time-stamping -- segment tree */
	m_segTree.initTimeStamping();

		/* run LK */
	m_segTree.setupCityOrder(m_tour);
	improved = 1;
    while( (t1 = m_lookbitQueue.deleteLookbit())>=0)
	{
        if( do_lk_search(t1, improved)>EPS)
		{
			improved++;
			m_segTree.setupCityOrder(m_tour);
		}
		else
			improved = 0;
	}

    trace("Quitting CLK::run\n");
}

	/*******************************************************************/
	/***************** LK algorithm ************************************/
	/* return gain. do_lk_search() implements step 2 ~ step 6 of LK paper. */
double CLK::do_lk_search(int t1, int improved)
{
		// t2 candidate list - defined two t2 about one t1 
	int             t2_cand_list[2];
	int             j1, j2, j3;
	
	_i = _k = 0;
	_G[0] = _Gstar = 0;
	t[1] = t1;
	
		/* get t2 candidate list */	
	m_tour->getNeighbor(t[1], &t2_cand_list[0], &t2_cand_list[1]);
	if( dist(t1, t2_cand_list[0])-EPS <= dist(t1, t2_cand_list[1]))
		SWAP(t2_cand_list[0], t2_cand_list[1], j3);
	
	for( j1=0; j1<2; j1++)
	{
		t[2] = t2_cand_list[j1];
		_i = 1;
		m_segTree.setupTree(1, t);
		
			/* get t3 candidate list */ 
		get_t3_cand_list();       // status : 0, seg(1)
		for( j2=0; j2<_num_t3_cand_list; j2++)
		{   // status : 0, seg(1)
			t[3] = _t3_cand_list[j2].t3;

				/* step 4. choose t4 to be a connected graph. normal t4 */
			t[4] = _t3_cand_list[j2].t4;
			_i = 2;
			make_two_change(improved);

			get_t5_cand_list();  /* get t5 cand. list with normal t4. */
			for( j3=0; j3<_num_t5_cand_list; j3++)
			{                    // status : 2change, seg(2)
				t[5] = _t5_cand_list[j3].t5;
				t[6] = _t5_cand_list[j3].t6;
				_i = 3;
				make_two_change(improved);

				if( _G[_i-1] > _Gstar+EPS)   /* further searching */
					search_deeper(improved); // status : 2-2change, seg(3)
				
				if( _Gstar > EPS) break;

					// when fail, reverse change to 2change, seg(2)
				m_tour->make2Change(t[1], t[6], t[5], t[4]);
				m_segTree.setupTree(2, t);
				_i = 2;
			}   // end of t5

			if( _Gstar > EPS) break; // status : 2change, seg(2)
			affirm(_k==0);
			
				// reverse change to 0, seg(1);
			m_tour->make2Change(t[1], t[4], t[3], t[2]);
			m_segTree.setupTree(1, t);
			_i = 1;
				
				/* try step 6.b - disconnected graph */
			t[4] = _t3_cand_list[j2].alter_t4; 
			_i = 2;
			get_alter_t5_cand_list();
			for( j3=0; j3<_num_t5_cand_list; j3++)
			{           // status : 0, seg(1);
				t[5] = _t5_cand_list[j3].t5;

					/* if t5 is between t2 and t3 */
				if( _t5_cand_list[j3].code == 1)
				{
					t[6] = _t5_cand_list[j3].t6;
					_i = 3;
					make_three_change(improved);
					search_deeper(improved); // status : 3change, seg(3)
					if( _Gstar>EPS) break;
						// reverse change to 0, seg(1);
					m_tour->make3Change(t[1], t[6], t[5],
										t[4], t[3], t[2]);
					m_segTree.setupTree(1, t);
					_i = 2;

					if( _t5_cand_list[j3].alter_t6>=0)
					{
						t[6] = _t5_cand_list[j3].alter_t6;
						_i = 3;
						make_three_change(improved); 
						search_deeper(improved);
						if( _Gstar>EPS) break;
							// reverse change to 0, seg(1)
						m_tour->make3Change(t[1], t[6], t[5],
											t[4], t[3], t[2]);
						m_segTree.setupTree(1, t);
						_i = 2;
					}
				}   // end of t6
				else /* if t5 lies between t1 and t4 */
				{
					t[6] = _t5_cand_list[j3].t6; _i = 3;
					t[7] = _t5_cand_list[j3].t7;
					t[8] = _t5_cand_list[j3].t8; _i = 4;
					make_four_change(improved);
					search_deeper(improved);
					if( _Gstar>EPS) break;
						// reverse change to 0, seg(1)
					m_tour->make4Change(t[1], t[8], t[7], t[6],
										t[5], t[4], t[3], t[2]);
					m_segTree.setupTree(1, t);
					_i = 2;
				}    // end of t6
			}   // end of alter t5
			if( _Gstar>EPS) break;
			_i = 1;
		}   // end of t3
		if( _Gstar>EPS) break;
	}   // end of t2
	if( _Gstar>EPS)
	{
		trace("t1=%d, j1=%d, j2=%d, alt j3=%d,", t1, j1, j2, j3);
		trace(" i=%d, k=%d, Gs=%d\n", _i, _k, _Gstar);
		reverse_change_to_best(improved);
	}
	else
	{
		trace("all backtracking failed - t1=%d\n", t1);
	}
	return _Gstar;
}
    /* search further. If not improved, reverse-change. */
void CLK::search_deeper(int improved)
{
    int         i, j, start_i;
    int         ci, ct, nt, nnt;
    int         best_nt, best_nnt;
	double Gix, gain, best_gain;

    start_i = _i;       // start_i is in no-change status.
    do
    {
            /* Get t[2i + 1], t[2i + 2] */
        ci = _i*2;  ct = t[ci];
        Gix = _G[_i-1] + dist(t[ci-1], ct);

        best_nt = INT_MIN;
		best_gain = -1e100;
        for( i=0; i<_nnn; i++)
        {
                // 1. Get a next t candidate of current t.
            nt = nni(ct, i);

                // 2. Check if yi(ct, nt) is valid.
                //    yi must be <not in T> and <not in {x1, x2, ..., xi-1}>.
            if( m_tour->isThereEdge(ct, nt)) continue;
                //    But In Johnson paper, he didn't use x-edge list.
                /*for( j=1; j<_i; j++)
                  if( m_tour->is_same_edge(ct, nt, t[j*2-1], t[j*2])) break;
                  if( j<_i) continue;*/
			
                // 3. Check if Gi > 0.
            if( Gix - dist(ct, nt) <=EPS) break;

                // 4. Get a next next t candidate.
			nnt = m_segTree.getPrev(nt);

                // 5. Check if x{i+1}(nt, nnt) are valid.
                //    1. Check if gain is better than best_gain.
                //    2. x{i+1}s must be <in T> and
                //       <not in {y1, y2, ..., yi-1}>.
			gain = dist(nt, nnt)-dist(ct, nt);
            if( gain <= best_gain+EPS) continue;

			for( j=1; j<_i; j++)
				if( m_tour->is_same_edge(nt, nnt, t[j*2], t[j*2+1]))
					break;
			if( j<_i) continue;

			best_nt = nt; best_nnt = nnt; best_gain = gain;
        }
        if( best_nt<0) break;

        t[ci+1] = best_nt;
        t[ci+2] = best_nnt;
        _i++;
		make_two_change(improved);
		
#if (LK_MORE_SPEED==1)
		if( _n>1000 && _i-_k>=10 && dist(t[ci+2], t[1]) > dist(t[ci], t[1])+EPS )
			break;
#endif
    } while( _i<=LK_DEPTH && _G[_i-1]>_Gstar+EPS);

    if( -EPS < _Gstar && _Gstar < EPS)
        reverse_change(start_i); // restore
}

int CLK::check_and_update(int improved)
{
    double gistar;
	
    gistar = dist(t[2*_i-1], t[2*_i]) - dist(t[2*_i], t[1]);
    if( _G[_i-1]+gistar > _Gstar+EPS)
    {
        _Gstar = _G[_i-1]+gistar;
        _k = _i;
        return 1;
    }
    return 0;
}

void CLK::reverse_change_to_best(int improved)
{
	int i;
	affirm(_k<=_i && _k>1);
	
	i = _i;
	while( i>_k)
	{
		i--;
		m_tour->make2Change(t[1], t[i*2+2], t[i*2+1], t[i*2]);
	}

	for( i=_k*2; i>0; i--) // set look-bit.
		m_lookbitQueue.addLookbit(t[i]);
}		

void CLK::reverse_change(int to)
{
	affirm(to<=_i && to>0);

	while( _i>to)
	{
		_i--;
		m_tour->make2Change(t[1], t[_i*2+2], t[_i*2+1], t[_i*2]);
	}
}

    /* Get (t3, t4) candidates */
void CLK::get_t3_cand_list()
{
    int         i, j, t3, t4, alter_t4;
	double x1_d;

    affirm( _i == 1);
    x1_d = dist(t[1], t[2]);

    _num_t3_cand_list = 0;
    for( i=0; i<_nnn; i++)
    {
            // 1. Get a t3 candidate.
        t3 = nni(t[2], i);

            // 2. Check if y1(t2, t3) is valid. y1 must be <not in T0>
		if( m_tour->isThereEdge(t[2], t3)) continue;

            // 3. Check if G1 > 0.
        if( x1_d - dist(t[2], t3) <=EPS) break;

            // 4. Get two t4 candidates.
		t4       = m_segTree.getPrev(t3);
		alter_t4 = m_segTree.getNext(t3);

            // 5. Check if x2s(t3, t4) are valid. x2s must be <in T0>.
            //    In this case, this condition is always satisified.
        _t3_cand_list[_num_t3_cand_list].t3 = t3;
		_t3_cand_list[_num_t3_cand_list].t4 = t4;
		_t3_cand_list[_num_t3_cand_list].alter_t4 = alter_t4;
		_t3_cand_list[_num_t3_cand_list].gain = dist(t3, t4) - dist(t[2], t3);
        _num_t3_cand_list++;
    }

        // Sort : maximize |x2| - |y1|.
	SORT(LK_T3_CAND, _t3_cand_list, _num_t3_cand_list,
		 _t3_cand_list[cur].gain > _t3_cand_list[key].gain+EPS);
}

    /* Get (t5, t6) candidates */
void CLK::get_t5_cand_list()
{
    int         i, j, t5, t6;
	double Gix;

    affirm(_i == 2);
    Gix = _G[1] + dist(t[3], t[4]);

    _num_t5_cand_list = 0;
    for( i=0; i<_nnn; i++)
    {
            // 1. Get a t5 candidate.
        t5 = nni(t[4], i);
		
            // 2. Check if y2(t4, t5) is valid.
            //    y2 must be <not in T1> and <not in {x1}>
        if( m_tour->isThereEdge(t[4], t5)) continue;
		if( m_tour->is_same_edge(t[4], t5, t[1], t[2])) continue;

            // 3. Check if G2 > 0.
        if( Gix - dist(t[4], t5) <=EPS) break;

            // 4. Get a t6 candidate.
		t6 = m_segTree.getPrev(t5);

            // 5. Check if x3(t5, t6) are valid.
            //    x3 must be <in T1> and <not in {y1}>.
		if( m_tour->is_same_edge(t5, t6, t[2], t[3])) continue;

        _t5_cand_list[_num_t5_cand_list].t5 = t5;
		_t5_cand_list[_num_t5_cand_list].t6 = t6;
		_t5_cand_list[_num_t5_cand_list].gain = dist(t5, t6) - dist(t[4], t5);
        _num_t5_cand_list++;
    }

        // sort : maximize |x3| - |y2|.
	SORT(LK_T5_CAND, _t5_cand_list, _num_t5_cand_list,
		 _t5_cand_list[cur].gain > _t5_cand_list[key].gain);
}

	/* Get (t5, t6) candidates when t4 is the abnormal point.*/
	/* We must solve only with m_tour and _city_order.
	   We must not use segment tree.
	   Careful for t[1] == t[4] case !!! */
void CLK::get_alter_t5_cand_list()
{
	int         i, j, t5, t6_1, t6_2, tmp;
	int         t6, t7, t8, t8_1, t8_2;
	int         best_t7, best_t8;
	double Gix, Gixx, gain, best_gain;
	
	affirm(_i == 2);
	Gix = _G[1] + dist(t[3], t[4]);
	
	_num_t5_cand_list = 0;
	for( i=0; i<_nnn; i++)
	{
			// 1. Get a t5 candidate.
		t5 = nni(t[4], i);
		
			// 2. Check if y2(t4, t5) is valid. 
			//    y2 must be <not in T1> and <not in {x1}>
        if( m_tour->isThereEdge(t[4], t5)) continue;
		if( m_tour->is_same_edge(t[4], t5, t[1], t[2])) continue;
		
			// 3. Check if G2 > 0.
		if( Gix - dist(t[4], t5) <=EPS) break;

			// 4. Get two t6 candidates.
		m_tour->getNeighbor(t5, &t6_1, &t6_2);
		
			// 5. Check if x3s(t5, t6) are valid. 
			//    x3s must be <in T1> and <not in {y1}>, <not in {x1, x2}.
		if( m_tour->is_same_edge(t5, t6_1, t[2], t[3]) ||
			m_tour->is_same_edge(t5, t6_1, t[1], t[2]) ||
			m_tour->is_same_edge(t5, t6_1, t[3], t[4])) t6_1 = -1;
		if( m_tour->is_same_edge(t5, t6_2, t[2], t[3]) ||
			m_tour->is_same_edge(t5, t6_2, t[1], t[2]) ||
			m_tour->is_same_edge(t5, t6_2, t[3], t[4])) t6_2 = -1;
		if( t6_1<0 && t6_2<0) continue;
		if( t6_1<0) SWAP(t6_1, t6_2, tmp);
		
			// 6. if t5 is between t2 and t3
		if( m_segTree.isBetween(t[2], t5, t[3]))
		{
			_t5_cand_list[_num_t5_cand_list].code     = 1;
				
				// 1. swap t6s such that maximizing |x3|.
			if( t6_2>=0 && dist(t5, t6_1) > dist(t5, t6_2)+EPS)
				SWAP(t6_1, t6_2, tmp);

			_t5_cand_list[_num_t5_cand_list].t6       = t6_1;
			_t5_cand_list[_num_t5_cand_list].alter_t6 = t6_2;
			_t5_cand_list[_num_t5_cand_list].gain = 
				dist(t5, t6_1) - dist(t[4], t5);
		}
		else
		{
			_t5_cand_list[_num_t5_cand_list].code = 0;

				// 1. Get valid t6.
			if( m_segTree.isBetween(t[4], t6_1, t5)) 
				t6 = t6_1;
			else if( t6_2>=0 && m_segTree.isBetween(t[4], t6_2, t5)) 
				t6 = t6_2;
			else continue;

			best_gain = -1e100;
			best_t7 = INT_MIN;
			Gixx = Gix - dist(t[4], t5) + dist(t5, t6);
			for( j=0; j<_nnn; j++)
			{
					// 1. Get a t7.
				t7 = nni(t6, i);
				
					// 2. Check if y3(t6, t7) is valid. 
					//    y3 must be <not in T1> and <not in {x1, y2}>
				if( m_tour->isThereEdge(t6, t7)) continue;
				if( m_tour->is_same_edge(t6, t7, t[1], t[2])) continue;
				if( m_tour->is_same_edge(t6, t7, t[4], t[5])) continue;
				
					// 3. Check if G3 > 0.
				if( Gixx - dist(t6, t7) <=EPS) 
					break;
				
					// 4. t7 must lie between t2 and t3.
				if( !m_segTree.isBetween(t[2], t7, t[3])) continue;

					// 5. Get two t8 candidates.
				m_tour->getNeighbor(t7, &t8_1, &t8_2);
				
					// 6. Check if x4s(t7, t8) are valid. 
					//    x4s must be <in T1> and <not in {y1, x1, x2}>.
				if( m_tour->is_same_edge(t7, t8_1, t[1], t[2]) ||
					m_tour->is_same_edge(t7, t8_1, t[3], t[4]))
					t8_1 = -1;
				if( m_tour->is_same_edge(t7, t8_2, t[1], t[2]) ||
					m_tour->is_same_edge(t7, t8_2, t[3], t[4]))
					t8_2 = -1;
					
				if( t8_1<0 && t8_2<0) continue;
				if( t8_1<0) SWAP(t8_1, t8_2, tmp);
				
					// 7. Get valid t8 that maximize |x4|.
				if( t8_2>=0 && dist(t7, t8_1) < dist(t7, t8_2)+EPS)
					t8 = t8_2;
				else
					t8 = t8_1;
				gain = dist(t7, t8) -  dist(t6, t7);

					// 8. compare current (t7, t8) with best (t7, t8)
				if( gain > best_gain+EPS)
				{
					best_t7 = t7; best_t8 = t8; best_gain = gain;
				}
			}
			if( best_t7 < 0) continue;
			
			_t5_cand_list[_num_t5_cand_list].t6 = t6;
			_t5_cand_list[_num_t5_cand_list].t7 = best_t7;
			_t5_cand_list[_num_t5_cand_list].t8 = best_t8;
			_t5_cand_list[_num_t5_cand_list].gain = best_gain +
				dist(t5, t6) -  dist(t[4], t5);
		}
		_t5_cand_list[_num_t5_cand_list].t5 = t5;
		_num_t5_cand_list++;
	}

		// sort : maximize gain.
	SORT(LK_T5_CAND, _t5_cand_list, _num_t5_cand_list,
		 _t5_cand_list[cur].gain > _t5_cand_list[key].gain+EPS);
}

	/*******************************************************************/
	/***************** Four Operations --- see Johnson paper ***********/
void CLK::make_two_change(int improved)
{   // Now, we know t[1] ~ t[2*i]. we know G[0] ~ G[i-2].
	int ci = (_i-1)*2;

		// Segment Tree
	m_segTree.do2Change(t[1], t[ci], t[ci+1], t[ci+2]);
	
		// C2EdgeTour
	m_tour->make2Change(t[1], t[ci], t[ci+1], t[ci+2]);

	_G[_i-1] = _G[_i-2] + dist(t[ci-1], t[ci]) - dist(t[ci], t[ci+1]);

	check_and_update(improved);
}

void CLK::make_three_change(int improved)
{   // Now, we know t[1] ~ t[6]. we know G[0].

		// Segment Tree
	m_segTree.do3Change(t[1], t[2], t[3], t[4], t[5], t[6]);
	
		// C2EdgeTour
	m_tour->make3Change(t[1], t[2], t[3], t[4], t[5], t[6]);

	_G[1] = dist(t[1], t[2]) - dist(t[2], t[3]);
	_G[2] = _G[1] + dist(t[3], t[4]) - dist(t[4], t[5]);

	check_and_update(improved);
}

void CLK::make_four_change(int improved)
{   // Now, we know t[1] ~ t[8]. we know G[0].

		// Segment Tree
	m_segTree.do4Change(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]);
	
		// C2EdgeTour
	m_tour->make4Change(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]);

	_G[1] = dist(t[1], t[2]) - dist(t[2], t[3]);
	_G[2] = _G[1] + dist(t[3], t[4]) - dist(t[4], t[5]);
	_G[3] = _G[2] + dist(t[5], t[6]) - dist(t[6], t[7]);
	
	check_and_update(improved);
}

	/*******************************************************************/
	/***************** Segment Tree class ******************************/
	/*******************************************************************/
CSegmentTree::CSegmentTree(int size)
{
	int i;

	_n                  = size;
	_city_order         = new int[_n];
	_time_stamp         = new int[_n];
	_loc_pointer        = new SEGMENT*[_n];
	
	_current_time = 1;  // initTimeStamping
	for( i=0; i<_n; i++) _time_stamp[i] = 0;
}

CSegmentTree::~CSegmentTree()
{
	delete[] _city_order;
	delete[] _time_stamp;
	delete[] _loc_pointer;
}

void CSegmentTree::initTimeStamping()
{
	/* 	
	int i;
	_current_time = 1;
	for( i=0; i<_n; i++) _time_stamp[i] = 0;
	*/
}

void CSegmentTree::setupCityOrder(C2EdgeTour* tour)
{
	affirm( tour && tour->isTour());
	int city, order;

	_tour = tour;
	_tour->findFirst(0);
	order = 0;
	while( (city = _tour->findNext()) >=0)
	{
		_city_order[city] = order++;
	}
	affirm( order==_n);
}
	
void CSegmentTree::setupTree(int to, int t[])
{
	affirm(to <= 2);
	
		// everey location pointers must be updated. change current time stamp.
	if( _current_time > INT_MAX - 100)
	{
		_current_time = 1;
	    for( int i=0; i<_n; i++) _time_stamp[i] = 0;
	}
	else
		_current_time++; 

	if( to<=2)      // make root
	{
		_seg_tree[0].first = 0; _seg_tree[0].last = _n-1;
		_seg_tree[0].left_child = _seg_tree[0].right_child = NULL;
		_seg_tree[0].prev = _seg_tree[0].next = &_seg_tree[0];
		if( (_city_order[t[1]]+1)%_n == _city_order[t[2]])
			_seg_tree[0].reversal = 0;
		else
			_seg_tree[0].reversal = 1;
		_seg_size = 1;
	}
	if( to==2)
	{
		affirm(_seg_size ==1);
		_tour->make2Change(t[1], t[4], t[3], t[2]);
		do2Change(t[1], t[2], t[3], t[4]);
		_tour->make2Change(t[1], t[2], t[3], t[4]);
	}
}	
	
CSegmentTree::SEGMENT* CSegmentTree::findSegment(int city)
{
	int      order;
	SEGMENT* seg = NULL;
	
	affirm(city>=0 && city < _n);


		// determine proper root to traverse the tree down.
	if( _time_stamp[city] == _current_time)
		seg = _loc_pointer[city];
	else
		seg = &_seg_tree[0]; // _seg_tree[0] is always real root of the tree.

		// traverse down
	order = _city_order[city];
	while( 1)
	{
		if( seg->first<=order && order<=seg->last)
		{
			if( _time_stamp[city] != _current_time)
			{
				_loc_pointer[city] = seg;
				_time_stamp[city]  = _current_time;
			}
			break;
		}
		else if( seg->first>order)
			seg = seg->left_child;
		else
		{
			affirm( seg->last<order);
			seg = seg->right_child;
		}
	}
	affirm( seg != NULL);
	return seg;
}

int CSegmentTree::getCityByOrder(int order, int referCity)
{
	int city;

	order = (order+_n)%_n;
	city = _tour->m_e1[referCity];
	if( city>=0 && _city_order[city] == order)
		return city;
	city = _tour->m_e2[referCity];
	if( city>=0 && _city_order[city] == order)
		return city;
	affirm(0); // invalid referCity or order;
	return -1;
}

void CSegmentTree::splitSegment(int c1, int c2)
{
	SEGMENT *s1, *s2;
	int     o1, o2, tmp;

	s1 = findSegment(c1);
	s2 = findSegment(c2);
	if( s1 != s2) return;
	
	SEGMENT *new1, *new2, seg;
	o1 = _city_order[c1]; o2 = _city_order[c2]; 
	if( o1 > o2)
	{
		SWAP(o1, o2, tmp); SWAP(c1, c2, tmp);
	}

	// rare case : when succ(c1)==c2, o1==_n-1 and o2=0
	if( o1==0 && o2==_n-1)
	{
		o1 = _n-2;
		c1 = getCityByOrder(o1, c2);
	}
	affirm( (o1+1)%_n == o2);

	seg = *s1;   // copy s1 to temp
	if( o1-seg.first >= seg.last-o2)
	{
		new1 = s1; new2 = &_seg_tree[_seg_size];
		new1->left_child = seg.left_child;
		new1->right_child = new2;
		new2->left_child = NULL;
		new2->right_child = seg.right_child;
	}
	else
	{
		new2 = s1; new1 = &_seg_tree[_seg_size];
		new1->left_child = seg.left_child;
		new1->right_child = NULL;
		new2->left_child = new1;
		new2->right_child = seg.right_child;
	}
	_seg_size++;
	
	new1->reversal = new2->reversal = seg.reversal;
	new1->first = seg.first; new1->last = o1;
	new2->first = o2;        new2->last = seg.last;
	
	if( seg.reversal == 0)
	{
		new1->next = new2;     new2->prev = new1;
		if( _seg_size == 2) // in case of first split, consider carefully
		{
			new1->prev = new2; new2->next = new1;
		}
		else
		{
			new1->prev = seg.prev; new2->next = seg.next; 
			seg.prev->next = new1;
			seg.next->prev = new2;
		}
	}
	else
	{
		new1->prev = new2;     new2->next = new1; 
		if( _seg_size == 2) // in case of first split, consider carefully
		{
			new1->next = new2; new2->prev = new1;
		}
		else
		{
			new1->next = seg.next; new2->prev = seg.prev;
			seg.prev->next = new2;
			seg.next->prev = new1;
		}
	}
	affirm( isCorrect());
}
		
void CSegmentTree::reverseTree()
{
	SEGMENT     *start, *end;
	
	if( _seg_size == 1)
		reverseSegment(&_seg_tree[0], &_seg_tree[0]);
	else if( _seg_size > 1)
	{
		end = &_seg_tree[0];
		start = _seg_tree[0].next;
		reverseSegment(start, end);
	}
}

void CSegmentTree::reverseSegment(SEGMENT* start, SEGMENT* end)
{
	SEGMENT *left, *right, *cur, *next, *tmp;
	
	if( start == end)
		start->reversal = !start->reversal;
	else if( start->prev == end)
	{
		affirm( end->next == start);
		
		cur = start;
		while( cur != end)
		{
			next          = cur->next;
			cur->reversal = !cur->reversal;
			SWAP(cur->next, cur->prev, tmp);
			cur           = next;
		}
		end->reversal = !end->reversal;
		SWAP(end->next, end->prev, tmp);
	}
	else
	{
		left = start->prev; right = end->next;
		left->next = end;   right->prev = start;

		cur = start;
		while( cur != end)
		{
			next          = cur->next;
			cur->reversal = !cur->reversal;
			SWAP(cur->next, cur->prev, tmp);
			cur           = next;
		}
		end->reversal = !end->reversal;
		SWAP(end->next, end->prev, tmp);
		start->next   = right;
		end->prev     = left;
	}
	affirm( isCorrect());
}

int CSegmentTree::isCorrect()
{
	SEGMENT *start, *current;
	int     num_segment;

	num_segment = 0;
	start = &_seg_tree[0];
	current = start->next; num_segment++;

	while( current != start)
	{
		num_segment++; current = current->next;
		affirm( num_segment <=(LK_DEPTH+2)*2);
	}
	affirm( num_segment == _seg_size);
	return (num_segment == _seg_size);
}

void CSegmentTree::do2Change(int t1, int t2, int t3, int t4)
{
	affirm( t2 == getNext(t1));
	affirm( t4 == getPrev(t3));
	SEGMENT *s1, *s2, *s3, *s4;
	
	splitSegment(t1, t2); splitSegment(t3, t4);
	
	s2 = findSegment(t2); s4 = findSegment(t4);
	affirm(s1 = findSegment(t1));
	affirm(s3 = findSegment(t3));
	affirm( s1->next == s2); affirm( s3->prev == s4);

	reverseSegment(s2, s4);
}

void CSegmentTree::do3Change(int t1, int t2, int t3, int t4, int t5, int t6)
{
	affirm( t2 == getNext(t1));
	affirm( t4 == getNext(t3));
	SEGMENT *s1, *s2, *s3, *s4, *s5, *s6;
	
	splitSegment(t1, t2); splitSegment(t3, t4);
	splitSegment(t5, t6);

	s2 = findSegment(t2); s3 = findSegment(t3);
	s5 = findSegment(t5); s6 = findSegment(t6);
	affirm( s1 = findSegment(t1));
	affirm( s4 = findSegment(t4));
	affirm( s1->next == s2); affirm( s3->next == s4);
	
	if( getNext(t5) == t6)
	{
		reverseSegment(s2, s3);
		reverseSegment(s3, s6); reverseSegment(s5, s2);
	}
	else
	{
		reverseSegment(s2, s6); reverseSegment(s5, s3);
	}
	affirm( findSegment(t3) == s3);
	affirm( findSegment(t6) == s6);
}

void CSegmentTree::do4Change(int t1, int t2, int t3, int t4, int t5, int t6,
					   int t7, int t8)
{
	affirm( t2 == getNext(t1));
	affirm( t4 == getNext(t3));
	affirm( t6 == getPrev(t5));
	SEGMENT *s2, *s3, *s4, *s6, *s7, *s8;
	
	splitSegment(t1, t2); splitSegment(t3, t4);
	splitSegment(t5, t6); splitSegment(t7, t8);

	s2 = findSegment(t2); s3 = findSegment(t3);
	s4 = findSegment(t4); s6 = findSegment(t6);
	s7 = findSegment(t7); s8 = findSegment(t8);
	
	if( getNext(t7) == t8)
	{
		reverseSegment(s2, s3); reverseSegment(s3, s8);
		reverseSegment(s7, s2); reverseSegment(s4, s6);
	}
	else
	{
		reverseSegment(s2, s8); reverseSegment(s7, s3);
		reverseSegment(s4, s6); 
	}
	affirm( findSegment(t3) == s3);
	affirm( findSegment(t8) == s8);
}
		
int CSegmentTree::getNext(int city)
{
	int c, a_order, b_order;
	SEGMENT *a, *b;

	a = findSegment(city);
	a_order = _city_order[city];
	
	if( a->reversal == 0)
	{
		if( a_order < a->last) b_order = a_order+1;
		else if( a_order == a->last)
		{
			b = a->next;
			b_order = (b->reversal==0) ? b->first : b->last;
		}
	}
	else
	{
		if( a_order > a->first) b_order = a_order-1;
		else if( a_order == a->first)
		{
			b = a->next;
			b_order = (b->reversal==0) ? b->first : b->last;
		}
	}
	c = getCityByOrder(b_order, city);
	affirm( _tour->isThereEdge(city, c));
	return c;
}
	
int CSegmentTree::getPrev(int city)
{
	int i, c, a_order, b_order;
	SEGMENT *a, *b;
	
	a = findSegment(city);
	a_order = _city_order[city];
	
	if( a->reversal == 0)
	{
		if( a_order > a->first) b_order = a_order-1;
		else if( a_order == a->first)
		{
			b = a->prev;
			b_order = (b->reversal==0) ? b->last : b->first;
		}
	}
	else
	{
		if( a_order < a->last) b_order = a_order+1;
		else if( a_order == a->last)
		{
			b = a->prev;
			b_order = (b->reversal==0) ? b->last : b->first;
		}
	}
	c = getCityByOrder(b_order, city);
	affirm( _tour->isThereEdge(city, c));
	return c;
}

int CSegmentTree::isBetween(int start, int middle, int end)
{
	if( middle == start || middle == end) return TRUE;
	if( start == end) return FALSE;
	
	int         starto, mido, endo, res;
	SEGMENT     *ss, *sm, *se, *cur, *next;
	
	starto = _city_order[start];
	mido   = _city_order[middle];
	endo   = _city_order[end];

	if( _seg_size == 1)
		ss = sm = se = &_seg_tree[0];
	else
	{
		ss = findSegment(start); 
		sm = findSegment(middle); se = findSegment(end);
	}
	
		// if ss == se, direction of path is important.
	if( ss == sm && ss == se)
	{
			// in this situation, we process 6 cases.
		if( starto < mido && mido < endo) res = TRUE;
		else if( starto < endo && endo < mido) res = FALSE;
		else if( mido < starto && starto < endo) res = FALSE;
		else if( mido < endo && endo < starto) res = TRUE;
		else if( endo < mido && mido < starto) res = FALSE;
		else res = TRUE;
		if( ss->reversal) res = !res;
	}
	else if( ss == se)
		res = (ss->reversal) ? (starto < endo) : (endo < starto);
	else if( ss == sm)
		res = (ss->reversal) ? (mido < starto) : (starto < mido);
	else if( se == sm)
		res = (sm->reversal) ? (endo < mido)   : (mido < endo);
	else   // now, ss != sm != se, thus calculate using segments.
	{
		cur = ss;
		while( (cur=cur->next) != se)
			if( cur == sm) break;
		res = (cur == sm);
	}
	return res;
}

	/*******************************************************************/
	/********************* CLookbitQueue Class *************************/
	/*******************************************************************/
CLookbitQueue::CLookbitQueue(int size)
{
	affirm( size > 0);
	m_look_bit_queue = new int[size];
	m_size           = size;
}

CLookbitQueue::~CLookbitQueue()
{
	delete[] m_look_bit_queue;
}

void CLookbitQueue::construct(C2EdgeTour* tour, 
						 C2EdgeTour* c1, C2EdgeTour* c2)
{
	affirm( tour->getSize() == m_size);
	
	int i, t;
	
	for( i=0; i<m_size; i++) m_look_bit_queue[i] = -1;
	m_look_head = m_look_tail = -1;
	if( c1==NULL)
	{
		tour->getNeighbor(m_size/2, &t, &i);
		for( i=0; i<m_size; i++)
			addLookbit((t+i)%m_size);
	}
	else
	{
		affirm( c2);
		for( i=0; i<m_size; i++)
		{
			if( !tour->isThereEdge(i, c1->m_e1[i]) ||
				!tour->isThereEdge(i, c1->m_e2[i]) ||
				!tour->isThereEdge(i, c2->m_e1[i]) ||
				!tour->isThereEdge(i, c2->m_e2[i]))
				addLookbit(i);
		}
	}
}

void CLookbitQueue::addLookbit(int c)
{
	affirm( 0<=c && c<m_size);
	
		/* There is already the city in look_queue */
	if( c == m_look_head || c == m_look_tail || m_look_bit_queue[c] >=0)
		return;

		/* case 1. when look_queue is empty. */
	if( m_look_tail<0)
		m_look_head = m_look_tail = c;
	else	/* case 2. when look_queue is not empty. */
	{
		m_look_bit_queue[m_look_tail] = c;
		m_look_tail = c;
	}
}

int CLookbitQueue::deleteLookbit()
{
	int res;
	if( m_look_head <0)        // When look_queue is empty.
		res = -1;
	else if( m_look_head==m_look_tail) // When there is one element in queue.
	{
		res = m_look_head;
		m_look_head = m_look_tail = -1;
	}
	else        // When there are two more elments in queue.
	{
		res = m_look_head;
		m_look_head = m_look_bit_queue[res];
		m_look_bit_queue[res] = -1;
	}
	return res;
}

	/*********************************************************************/
	/********************* Test Driver ******* ***************************/
#ifdef LK_TEST
#include "cpputil.h"

void setupSS();
void goSS(int step);
int main(int argc, char* argv[])
{
	CLK*            lk;
	C2EdgeTour*     tour;
	int *OptPath;
	int n;

	ReadTspFile((char *)"cycle.in"); ConstructNN(20);

	lk   = new CLK(gNumCity, gNumNN);
	tour = new C2EdgeTour(gNumCity);
	
	srandom(time(0));

	n = gNumCity;


	OptPath = new int[n];
	tour->makeRandomTour();
	lk->run(tour);
	double tourcost = tour->evaluate();
	tour->convertToOrder(OptPath, n);

	FILE *fout = fopen("cycle.out","w");
	//fprintf(fout,"%lf\n",tourcost);

	/*
	 * old:
	 * for(int i=1;i<n;i++) fprintf(fout,"%d ",OptPath[i]); fprintf(fout,"\n");
	 */
	for(int i=0;i<n;i++) fprintf(fout,"%d ",OptPath[i]+1); fprintf(fout,"\n");

	fclose(fout);

	delete lk; delete tour; delete OptPath;
	return 0;
}

int (*gOrigDist)[20];
int gMinDist;
void setupSS()
{
	int i, j;

	if( gOrigDist) delete[] gOrigDist;
	gOrigDist = new int[gNumCity][20];
	
	gMinDist = 200000000;
	for( i=0; i<gNumCity; i++)
		for( j=0; j<2; j++)
		{
			gOrigDist[i][j] = nnd(i, j);
			if( nnd(i, j) < gMinDist) gMinDist = nnd(i, j);
		}
	// space smoothing 
	for( i=0; i<gNumCity; i++)
		for( j=0; j<2; j++)
		{
			nnd(i, j) = gMinDist; // It should work!!
		}
}
			
void goSS(int reverse)
{
	int step = 1 - reverse;
	int i, j;

	for( i=0; i<gNumCity; i++) // recover 'step' columns
	{
		nnd(i, step) = gOrigDist[i][step];
	}
}
	

#endif








