#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "macro.h"
#include "tsplib_io.h"
#include "tour.h"

	////////////////////////////// DEBUG /////////////////////////////////
#define affirm(x)
#define trace

//#define affirm(x)    assert(x)
	//#define trace        printf

	/**********************************************************************/
	/****************************  CTour Implementation *****************/
CTour::CTour()
{
	m_length = -1;
	m_size   = 0;
}

int CTour::getHammingDist(CTour* tour)
{
	int v1, v2, count;

	enumEdgeFirst(0);
	count = 0;
	while( enumEdgeNext(&v1, &v2))
		if( tour->isThereEdge(v1, v2)) count++;
	affirm( getSize() - count >= 0);
	return (getSize() - count);
}

double CTour::evaluate()
{
	int v1, v2;
	m_length = 0;

	enumEdgeFirst(0);
	while( enumEdgeNext(&v1, &v2))
		m_length += dist(v1, v2);
	return m_length;
}
			
void CTour::enumEdgeFirst(int start)
{
	findFirst(start);
	m_edge_start = m_edge_v1 = findNext();
}

int CTour::enumEdgeNext(int* v1, int* v2)
{
	if( m_edge_start < 0) return FALSE;
	*v1 = m_edge_v1;
	m_edge_v1 = findNext();
	*v2 = m_edge_v1;
	if( m_edge_v1 < 0)
	{
		*v2 = m_edge_start; m_edge_start = -1;
	}
	affirm( isThereEdge(*v1, *v2));
	return TRUE;
}

void CTour::convertToOrder(int* dest, int size)
{
	affirm( size >= gNumCity);
	int     city, order;

	findFirst(0);
	order = 0;
	while( (city=findNext()) >= 0)
		dest[order++] = city;
}

void CTour::makeRandomOrderTour(int* dest, int size)
{
	int     i, j, k, t;
	for( i=0; i<size; i++) dest[i] = i;
	
	for( i=0; i<size/2; i++)
	{
		j = random()%size; k = random()%size;
		SWAP(dest[j], dest[k], t);
	}
}

	/*********************************************************************/
	/************************ C2EdgeTour Implementation ******************/
C2EdgeTour::C2EdgeTour(int size)
{
	create(size);
}

C2EdgeTour::C2EdgeTour()
{
	m_size = 0; m_e1 = m_e2 = NULL;
}

C2EdgeTour::~C2EdgeTour()
{
	delete[] m_e1;
	delete[] m_e2;
}

void C2EdgeTour::create(int size)
{
	affirm(size>0);
	
	int i;
	
	m_size = size;
	m_e1 = new int[m_size];
	m_e2 = new int[m_size];
	m_isMult = FALSE;
	for( i=0; i<m_size; i++)
		m_e1[i] = m_e2[i] = -1;
}

void C2EdgeTour::makeRandomTour()
{
	int i, *order_array;
	
	order_array = new int[m_size];
	makeRandomOrderTour(order_array, m_size);
	for( i=0; i<m_size; i++)
		m_e1[i] = m_e2[i] = -1;
	
	for( i=0; i<m_size-1; i++)
		addEdge(order_array[i], order_array[i+1]);
	addEdge(order_array[0], order_array[m_size-1]);
	delete[] order_array;
	affirm(isTour());
}

void C2EdgeTour::findFirst(int start)
{
	m_prev = m_end = m_e2[start]; m_cur = start;
}

int C2EdgeTour::findNext()
{
	int temp;
	if( m_cur < 0) return -1;
	temp = m_cur;
	m_cur = ( m_e1[m_cur] == m_prev) ? m_e2[m_cur] : m_e1[m_cur];
	m_prev = temp;
	if( m_prev == m_end) m_cur = -1;
	return m_prev;
}
	
int C2EdgeTour::isThereEdge(int v1, int v2)
{
	affirm( ((m_e1[v1] == v2)||(m_e2[v1] == v2)) ==
			   ((m_e1[v2] == v1)||(m_e2[v2] == v1)));
	return ((m_e1[v1] == v2) || (m_e2[v1] == v2));
}


void C2EdgeTour::addEdge(int v1, int v2)
{
	affirm( m_isMult || !isThereEdge(v1, v2));
	affirm(v1 != v2);
	affirm(m_e1[v1]<0 || m_e2[v1]<0);
	affirm(m_e1[v2]<0 || m_e2[v2]<0);
		//trace("In addEdge, v1=%d, v2=%d\n", v1, v2);
	
		/* add v2 in v1 side */
	if( m_e1[v1]<0) m_e1[v1] = v2;
	else            m_e2[v1] = v2;
	
		/* add v1 in v2 side */
	if( m_e1[v2]<0) m_e1[v2] = v1;
	else            m_e2[v2] = v1;
}

void C2EdgeTour::addEdgeSafely(int v1, int v2)
{
	if( !isThereEdge(v1, v2)) addEdge(v1, v2);
}

void C2EdgeTour::deleteEdge(int v1, int v2)
{
	int *t;
	affirm( isThereEdge(v1, v2));
	
		/* delete v2 in v1 side */
	affirm(v1 != v2);
	if( m_e1[v1]==v2) m_e1[v1] = -1;
	else              m_e2[v1] = -1;

		/* delete v1 in v2 side */
	if( m_e1[v2]==v1) m_e1[v2] = -1;
	else              m_e2[v2] = -1;
}

int C2EdgeTour::getNumEdge()
{
	int i, n=0;

	for( i=0; i<m_size; i++)
	{
		if( m_e1[i] >=0) n++;
		if( m_e2[i] >=0) n++;
	}
	affirm( n%2 == 0);
	return (n/2);
}

int C2EdgeTour::isTour()
{
	int i, n;
	
	for( i=0; i<m_size; i++) if( getEdgeSize(i) <= 1) return FALSE;

	findFirst(0);
	n = 0;
	while( findNext()>=0)
	{
		if( n > m_size) return FALSE;
		n++;
	}
	if( n < m_size) return FALSE;
	return (getNumEdge() == m_size);
}

int C2EdgeTour::isDisjointCycle()
{
	int i, n;
	for( i=0; i<m_size; i++)
	{
		affirm( getEdgeSize(i)==2);
	}

	findFirst(0);
	n =0;
	while( findNext()>=0)
		n++;
	if( n == m_size) return FALSE;
	else return TRUE;
}

void C2EdgeTour::make2Change(int t1, int t2, int t3, int t4)
{
	deleteEdge(t1, t2); deleteEdge(t3, t4);
	addEdge(t1, t4);    addEdge(t2, t3);
}

void C2EdgeTour::make3Change(int t1, int t2, int t3, int t4,
							 int t5, int t6)
{
	deleteEdge(t1, t2); deleteEdge(t3, t4); deleteEdge(t5, t6);
	addEdge(t1, t6);    addEdge(t2, t3);    addEdge(t4, t5);
}

void C2EdgeTour::make4Change(int t1, int t2, int t3, int t4,
							 int t5, int t6, int t7, int t8)
{
	deleteEdge(t1, t2); deleteEdge(t3, t4); 
	deleteEdge(t5, t6); deleteEdge(t7, t8);
	addEdge(t1, t8);    addEdge(t2, t3);
    addEdge(t4, t5);    addEdge(t6, t7);
}
	
int C2EdgeTour::isBetween(int prev, int cur, int city, int end)
{
	if( city == cur || city == end) return TRUE;

	int next;
	while( (next = getNext(prev, cur)) >= 0)
	{
		if( next == end) return FALSE;
		if( next == city) return TRUE;
		prev = cur;
		cur  = next;
	}
	return FALSE;
}

int C2EdgeTour::convertFromOrder(int* src, int size)
{
	affirm(size == m_size);
	int i;
	
	for( i=0; i<m_size; i++) m_e1[i] = m_e2[i] = -1;
	for( i=0; i<m_size; i++)
		addEdge(src[i], src[(i+1)%m_size]);
	affirm(isTour());
}

const C2EdgeTour& C2EdgeTour::operator=(const C2EdgeTour& src)
{
	affirm(src.m_size == m_size);

	memcpy(m_e1, src.m_e1, sizeof(int)*m_size);
	memcpy(m_e2, src.m_e2, sizeof(int)*m_size);
	m_length = src.m_length;
	return *this;
}

void C2EdgeTour::constructTabuEdge(C2EdgeTour* tabu)
{
	int     i, ep;
	
	tabu->deleteAllEdge();
	for( i=0; i<m_size; i++)
	{
		ep = getEndPoint(i);
		if( i < ep) tabu->addEdge(i, ep);
	}
}

int C2EdgeTour::canAddEdge(C2EdgeTour* tabu, int v1, int v2)
{
	if( tabu->isThereEdge(v1, v2) || isThereEdge(v1, v2) ||
		getEdgeSize(v1) == 2 || getEdgeSize(v2) == 2)
		return FALSE;
	else
		return TRUE;
}

void C2EdgeTour::addEdgeWithTabu(C2EdgeTour* tabu, int v1, int v2)
{
	affirm( tabu->getEdgeSize(v1)<2 && tabu->getEdgeSize(v2)<2);
	int     n1, n2, ep1, ep2;
	
	addEdge(v1, v2);
	
	tabu->getNeighbor(v1, &n1, &n2);
	ep1 = (n1>=0) ? n1 : n2;
	if( ep1 >=0) tabu->deleteEdge(v1, ep1);
	else         ep1 = v1;
	
	tabu->getNeighbor(v2, &n1, &n2);
	ep2 = (n1>=0) ? n1 : n2;
	if( ep2 >=0) tabu->deleteEdge(v2, ep2);
	else         ep2 = v2;
	
	tabu->addEdge(ep1, ep2);
}

void C2EdgeTour::connectRandom(C2EdgeTour* tabu) 
{
	int         i, j, k, nv, *av, size, rne, tt;
	
	size = getSize();
	av = new int[size];
	
		/* save the cities that there are in the end of a path */
	nv = 0;
	for( i=0; i<size; i++)
		if( getEdgeSize(i) <= 1)
			av[nv++] = i;

	for( i=0; i<nv/2; i++)
	{
		j = random()%nv; k = random()%nv;
		SWAP(av[j], av[k], tt);
	}
	rne = size - getNumEdge(); // we must add the rne number of edges.
	while( rne>1)
	{
		for( i=0; i<nv; i++)
			if( av[i]>=0) break;
		j = i;
		for( i++; i<nv; i++)
			if( av[i]>=0 && !tabu->isThereEdge(av[j], av[i])) break;
		k = i;
		affirm( av[j]>=0 && av[k]>=0);
		affirm( getEdgeSize(av[j])<2 && getEdgeSize(av[k])<2);
		
		addEdgeWithTabu(tabu, av[j], av[k]);
		if( getEdgeSize(av[j]) == 2) av[j] = -1;
		if( getEdgeSize(av[k]) == 2) av[k] = -1;
		
		rne--;
	}
		/* necessary step. 
		   The use of tabu edge make the tour a Hamiltonian PATH after all,
		   less the tour is already a complete tour. (before this 
		   function call)  */
	affirm( getNumEdge() >= getSize()-1);
	if( getNumEdge() == getSize()-1)
	{
		for( i=0; i<nv; i++)
			if( av[i] >= 0) break;
		j = i;
		for( i++; i<nv; i++)
			if( av[i] >= 0) break;
		k = i;
		addEdge(av[j], av[k]);
	}
	affirm( isTour());
	delete[] av;
}

int C2EdgeTour::getEndPoint(int v)
{
	int         es, ep, prev_ep, t;

	es = getEdgeSize(v);

	if( es == 0 || es == 2) return -1;

	ep = ( m_e1[v] >= 0) ? m_e1[v] : m_e2[v];
	prev_ep = v;

	while( getEdgeSize(ep) == 2)
	{
		t  = ep;
		ep = (m_e1[ep]==prev_ep) ? m_e2[ep] : m_e1[ep];
		prev_ep = t;
		if( ep == v) return -1;
	}
	affirm( getEdgeSize(ep) == 1);
	return ep;
}


static int compareMin(const void* a1, const void* a2)
{
	IntTriple* b1 = (IntTriple*)a1;
	IntTriple* b2 = (IntTriple*)a2;
	
	if( b1->i3 < b2->i3) return (-1);
	if( b1->i3 > b2->i3) return 1;
	return 0;
}


















































































