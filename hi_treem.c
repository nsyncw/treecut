/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>

#include "types_treem.h"  /* type definitions */
#include "parser_treem.c" /* parser */
#include "timer.c"        /* timing routine */

/*
#define GLOB_UPDT_FREQ 0.5
*/
#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2

#define MAX_PASSES 3000
#define MAX_LONG LONG_MAX
NodePropArr allResults[MAX_PASSES];
long N, M;
nodeP *nodes = NULL;

cType *gpfa = NULL;
cType *gpdep = NULL;
long *gpcv = NULL;
short *gps = NULL;
cType *gpaccup = NULL;
long *gpaccmcv = NULL;
cType *gpaccposmcv = NULL;
cType *gpaccjointid = NULL;
long *gpaccjointmcv = NULL;

int mode;

//----------------FUNCTION DEFINITION:

void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}

cType *randNums = NULL;
cType randNumIdx = 0;
cType mrand()
{
  return randNums[randNumIdx++];
}

void initrand()
{
  //(1) init rand

  srand((int)(timer() * 1000));
  if (randNums == NULL)
  {
    randNums = (cType *)walloc(2 * M, sizeof(cType));
  }

  randNumIdx = 0;
  for (int i = 0; i < 2 * M; i++)
  {
    randNums[i] = (cType)rand();
    // if(i<100){
    //   printf("rand %d %ld\n",i,randNums[i]);
    // }
  }
}

/////////////////////////////quicksort
//划分数组的函数
int split(sType* idx, edgeP* edges, int low, int high)
{
    sType t;
    int i = low;    //i指向比较元素的期望位置
    int x = edges[idx[i]].tmp;    //将该数组第一个元素设置为比较元素
    //从数组的第二个元素起开始遍历，若找到的元素大于比较元素，则跳过
    for(int j = low+1;j<=high;j++){
        //若找到了小于比较元素的数，则将其与前面较大的数进行交换
        if (edges[idx[j]].tmp < x)
        {
            i++;
            t = idx[i]; idx[i] = idx[j]; idx[j] = t;
        }
    }

    t = idx[i]; idx[i] = idx[low]; idx[low] = t;    //将比较元素交换到期望位置
    return i;
}

//快速排序
void quicksort(sType *idx, edgeP * edges , int low, int high)
{
    if (low < high)
    {
        int i = split(idx, edges, low, high);    //划分数组并获得比较元素位置
        quicksort(idx, edges,  low, i - 1);    //对比较元素左边进行排序
        quicksort(idx, edges,  i + 1, high);    //对比较元素右边进行排序
    }
}


/////////////////////////////////////heap sort
void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    // assert(idx[0] != idx[3]);
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];
        // printf("start %d <= i %d \n",start,i);

        // if(start == 3 || start == 2 || i==2 || i==3){
        //     assert(idx[0] != idx[3]);
        // }

        start = i;  
        i = 2*start+1;  
    }  

    // printf("fin start %d <= tmp %d \n",start,tempIdx);
    idx[start] = tempIdx;  

    // for(int k=0; k<10; k++){
    //     printf("%d, ",idx[k]);
    // }
    // printf("\n");

    // assert(idx[2] != idx[3]);
    // assert(idx[0] != idx[3]);
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        // printf("swap 0 with %d \n",i);
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        // if( i==2 || i==3){
        //     assert(idx[0] != idx[3]);
        //     assert(idx[2] != idx[3]);
        // }

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

    // assert(idx[2] != idx[3]);
}  
  
// void percdown(sType *idx, edgeP * edges , int i, int n)
// {
// 	int child;        //孩子结点
// 	int temp;		  //临时变量
// 	/*对堆进行整理*/                                       
// 	for(temp = i; i * 2 + 1 < n; i = child)
// 	{
// 		child = i * 2 + 1; 								//当前结点的左孩子结点
// 		if(child  != n - 1 && edges[idx[child]].tmp < edges[idx[child+1]].tmp)  //比较左孩子和右孩子谁大
// 			child++;
// 		if(edges[idx[temp]].tmp < edges[idx[child]].tmp)
// 			idx[i] = idx[child];
// 		else
// 			break;
// 	}

// 	idx[i] = idx[temp]; 
// }


// void heapsort(sType *idx, edgeP * edges , int n)   //n为数组长度
// {
// 	int i;
// 	/*初始化堆*/
//     sType tmp;
// 	for(i = n / 2; i >= 0; i--)
// 		percdown(idx,edges, i, n);
// 	/*排序*/	
// 	for(i = n - 1; i >=0 ; i--)
// 	{
// 		tmp = idx[0];
//         idx[0] = idx[i];
//         idx[i] = tmp;           //将堆顶元素沉入堆底，并且将堆长度减1
// 		percdown(idx,edges, 0, i); 
// 	}
// }

//desending
void deOrderEdgeByRandomCap(nodeP *np)
{

  cType cnt = np->nIdx;
  //assert(cnt < 1000);

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    pedges[i].tmp = -1*mrand() % pedges[i].cap;
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
//   quicksort(idxs,pedges,0,cnt-1);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}

//ascending
void aOrderEdgeByWeight(nodeP *np)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    pedges[i].tmp = mrand() % pedges[i].w; 
  }

  for (int i = 1; i < cnt; i++)
  {
    for (int j = i; j > 0; j--)
    {
      if (pedges[idxs[j]].tmp < pedges[idxs[j - 1]].tmp) 
      {
        sType tmp = idxs[j];
        idxs[j] = idxs[j - 1];
        idxs[j - 1] = tmp;
      }
    }
  }
}

void aOrderEdgeByAvgCV(nodeP *np)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand() % acv; 
    }
  }

    // quicksort(idxs,pedges,0,cnt-1);
    HeapSort(idxs,pedges,cnt);
}

// void deOrderEdgeByRandomCap(edgeP** pedges){
//   edgeP* edges = *pedges;
//   edgeP* eh = edges;
//   edgeP* nextEh = edges->next;
//   //pop order
//   while(eh != NULL){
//     // printf("cur endNode %d\n",eh->endNode);
//     eh->tmp = rand()%(eh->cap*10);
//     for(edgeP* nnh = eh; nnh->prev != NULL; ){

//       edgeP* pnh = nnh->prev;
//       if(nnh->tmp > pnh->tmp ){

//         edgeP* tmpp = pnh->prev;
//         edgeP* tmpn = nnh->next;

//         pnh->next = tmpn;
//         nnh->prev = tmpp;
//         pnh->prev = nnh;
//         nnh->next = pnh;

//         if(tmpp == NULL){
//           *pedges = nnh;
//         }
//         else{
//           tmpp->next = nnh;
//         }

//         if(tmpn != NULL){
//           tmpn->prev = pnh;
//         }

//       }
//       else{
//         nnh = nnh->prev;
//       }
//     }

//     eh = nextEh;
//     if(eh != NULL){
//       nextEh = eh->next;
//     }
//   }
// }

//ascending

#define min(s, t) (s) < (t) ? (s) : (t)
#define max(s, t) (s) > (t) ? (s) : (t)

//build acc for one time
//upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
cType SPAN_LEN=0;
void buildAcc(cType curN, cType upid, long mcv, cType lastDepMCV,cType lastJointNodeId, long lastJointMCV)
{
  assert(curN <= N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  edgeP *pedges = (nodes + curN)->edges;
  long curCV = gpcv[curN];
  cType curDep = gpdep[curN];

  assert(gpdep[curN] == 0 || curCV >0);
//logic for span
  mcv = min(mcv, curCV);
  if(mcv == curCV){
    lastDepMCV = curDep;
  }

  gpaccmcv[curN] = mcv;
  gpaccup[curN] = upid;
  gpaccposmcv[curN] = curDep - lastDepMCV;

  if (curDep % SPAN_LEN == 0)
  {
    upid = curN;
    mcv = MAX_LONG;
    lastDepMCV = 0; //doesn't matter, will be udpated in the subsequent call
  }

//-logic for joint node
  cType childCnt = gpaccjointid[curN];
  lastJointMCV = min(lastJointMCV,curCV);
  gpaccjointmcv[curN] = lastJointMCV;
  gpaccjointid[curN] = lastJointNodeId;

  assert(gpdep[curN] == 0 || gpaccjointmcv[curN] > 0.1);
  
  //gpaccjointid[curN] is updated by markcut() to contain the number of traversing children
  if(childCnt > 1 ){
    lastJointNodeId = curN;
    lastJointMCV = MAX_LONG;
  }
  else{
    //do nothing
  }


  while (cnt > 0)
  {
    if (gpfa[pedges->endNode] == curN)
    {
      buildAcc(pedges->endNode, upid, mcv,lastDepMCV,lastJointNodeId,lastJointMCV);
    }
    pedges++;
    cnt--;
  }

  // assert(childCnt == 0);
}

//VER 1
// void buildAcc(cType curN, cType upid, long mcv)
// {
//   assert(curN <= N && curN >= 1);
//   int cnt = (nodes + curN)->nIdx;
//   edgeP *pedges = (nodes + curN)->edges;

//   mcv = min(mcv, gpcv[curN]);
//   gpaccmcv[curN] = mcv;
//   gpaccup[curN] = upid;

//   if (gpdep[curN] % SPAN_LEN == 0)
//   {
//     upid = curN;
//     mcv = MAX_LONG;
//   }

//   while (cnt > 0)
//   {
//     if (gpfa[pedges->endNode] == curN)
//     {
//       buildAcc(pedges->endNode, upid, mcv);
//     }
//     pedges++;
//     cnt--;
//   }
// }

// //w can be propagate for two degree node
// //logic is the neighboring two degree edge 's w is the accumulation of all this edges
// //if it is -<, edge in < should has w for cnt+1 cut set
// //w propagate should be right
// void reupdateWeight(){

// }

void preMarkCut(cType curN)
{
  assert(curN <= N && curN >= 1);
  nodeP * np = nodes+curN;
  int cnt = np->nIdx;
  edgeP *pedges = np->edges;
  cType tmp = 0;
  while (cnt > 0)
  {
    tmp += pedges->cap;
    pedges++;
    cnt--;
  }

  np->totalCap = tmp;

}


//mark for one time;
void markCut(cType curN)
{
  assert(curN <= N && curN >= 1);

  assert((nodes + curN)->nIdx > 0);

  short *curS = gps + curN;

  assert(*curS == 0);

  long *curCV = gpcv + curN;
  cType *curDep = gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedEdges;

  // //DEBUG
  // for(int i=0; i<cnt; i++){
  //   if(idxs[i] >= cnt){
  //     assert(idxs[i] < cnt);
  //   }
  // }

  if (mode == 1)
  {
    // printf("cur node %d\n",curN);
    deOrderEdgeByRandomCap(np); //
  }
  else if (mode == 2)
  {
    // aOrderEdgeByWeight(np); //
    aOrderEdgeByAvgCV(np);
  }

  //   //DEBUG
  // for(int i=0; i<cnt; i++){
  //   assert(idxs[i] < cnt);
  // }
  cType uplinkcap = 0;
  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = gps[zn];

    assert(!(gpfa[zn] == curN && zs == 2));


    if (zs == 1)
    {
      cap = eh->cap;
      *curCV += cap;
      gpcv[zn] -= cap;

      //not incluing uplink, since the value 1 of uplink is useless, each node has it
      if (mode == 1 && gpdep[zn] != *curDep - 1)
      {
        uplinkcap += eh->cap;
      }

      // if (mode == 1 && gpdep[zn] != *curDep - 1)
      // {
      //   eh->w += 1;
      //   eh->rev->w += 1;
      // }
    }
    else if (zs == 0)
    {
      gpfa[zn] = curN;
      gpdep[zn] = *curDep + 1;
      gpaccjointid[curN] ++;
      markCut(zn);
      assert(gpdep[zn] == gpdep[curN] + 1);
      *curCV += gpcv[zn];
    }
    else
    {
      //bypass, no need to handle
    }
  }

  short allset = 0;
  if((np->totalCap - uplinkcap) <= uplinkcap){
    allset = 1;
  }

  if(mode == 1){
//update ver 2 w,according to
    for (int ni = 0; ni < cnt; ni++)
    {
      // nodeP* znp = nodes+eh->endNode;

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      // nodeP *znp = nodes+zn;

      assert(zn != 0);
      assert(zn != curN);
      short zs = gps[zn];
    //   assert(zs != 2);
      //propgate weight to zn's edges
      // if(zs == 1 && gpdep[zn] != *curDep - 1){
      //   if(znp->totalCap - eh->cap <= eh->cap){
      //     edgeP* zep = znp->edges;
      //     int tmpk = znp->nIdx;
      //     while(tmpk >0){
      //       cType weight = zep-> w;

      //       if(zep->avgCV == 0){
      //         zep->avgCV = MAX_LONG;
      //       }
      //       zep->avgCV = min(zep->avgCV, *curCV);//((eh->avgCV) * weight + *curCV)/(weight+1);
      //       zep->w = weight+1;
            
      //       edgeP *rzep = zep->rev;

      //       if(rzep->avgCV == 0){
      //         rzep->avgCV = MAX_LONG;
      //       }

      //       weight = rzep-> w;
      //       rzep->avgCV = min(rzep->avgCV, *curCV);//((reh->avgCV) * weight + *curCV)/(weight+1);
      //       rzep->w = weight+1;
            
      //       zep ++;
      //       tmpk--;
      //     }
      //   }

      // }

      //progate weight to curN's edges
      if (allset ==1 || (allset == 0 && zs == 1 && gpdep[zn] != *curDep - 1))
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);//((eh->avgCV) * weight + *curCV)/(weight+1);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);//((reh->avgCV) * weight + *curCV)/(weight+1);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(gpdep[curN] == 0 || *curCV >0);

  *curS = 2;

  // assert(allResults[0].pfa[0] == 0);
}

// void mcvassert(cType *pfa, long *pCV, cType *paccup, long *paccmcv, cType n)
// {
//   long mcv = MAX_LONG;
//   cType upn = paccup[n];
//   long omcv = paccmcv[n];
//   while (n != upn)
//   {
//     mcv = min(mcv, pCV[n]);
//     n = pfa[n];
//   }

//   assert(mcv == omcv);
// }

//mincandi is the current smallest value by other passes, used to accelerate
long solveMaxFlowAccVER4(long minCandi, NodePropArr np, cType s, cType t)
{
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;
  cType *paccup = np.pacc_upid;
  long *paccmcv = np.pacc_upmincv;
  cType *paccposmcv = np.pacc_pos_upmincv;
  cType *jup = np.pacc_jointid;
  long *jmcv = np.pacc_jointmcv;

  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  cType depT = pDep[t];
  assert(pDep[s] >= depT);

  long mcv = MAX_LONG;
  cType ups = paccup[s];

  while (pDep[ups] > depT)
  {
    assert(pDep[ups] % SPAN_LEN == 0);
    mcv = min(mcv, paccmcv[s]);
    s = ups;
    ups = paccup[ups];
  } 
//   assert(mcv >100.0);
  assert(pDep[ups] <=depT && pDep[s] >= depT);

  cType upt = t;

  if (pDep[t] % SPAN_LEN != 0)
  {
    upt = paccup[t];
  }

  assert(pDep[ups] == pDep[upt]);

  while (ups != upt)
  {
    mcv = min(mcv, min(paccmcv[s], paccmcv[t]));

    s = ups;
    ups = paccup[ups];

    t = upt;
    upt = paccup[upt];

    assert(pDep[s] % SPAN_LEN == 0);
    assert(pDep[t] == pDep[s]);
  }        
  
  assert(ups == upt);
  if(s == t){
    return mcv;
  }
  
  cType min_bound2 = min(paccmcv[s],paccmcv[t]);

  if(min_bound2 >= mcv || min_bound2 >= minCandi){
    //no need to search
    return mcv;
  }

  if(min_bound2 == paccmcv[s] && pDep[t]+paccposmcv[s] < pDep[s]){
    return mcv;
  }


  ///////////////////////check inside one SPAN
  //(1) we need to check whether s and t in the same line, i.e., t is the ancestor of s  
  //the only way is to apprach the depth of t and check

  if(pDep[s] != pDep[t]){

    if(pDep[s] < pDep[t]){
      cType tmp = s;
      s = t;
      t = tmp;
    }  

    cType jups = jup[s];
    cType depT = pDep[t];
    while(pDep[jups] > depT){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }

      s = jups;
      jups = jup[s];

    }

    if(pDep[jups] == depT){
      if(jups == t){
        return min(mcv,jmcv[s]);
      }
      else{
        //s and t in two lines
        assert(pDep[jups] == depT && pDep[s] > depT);
        goto STEP_CHECK_IN_TWOLINES;
      }
    }
    else{
      assert(pDep[jups] < depT && pDep[s] > depT);
      while(pDep[s] > depT){
        mcv = min(mcv,pCV[s]);
        if(mcv == min_bound2){
          return mcv;
        }
        s = pFa[s];        
      }

      if(s == t){
        return mcv;
      }

      assert(s!=t && pDep[s] == depT);
      //s and t in two lines.
      goto STEP_CHECK_IN_TWOLINES;

    }
  }
  else{
      //pDep[s] == pDep[t];
      if(s == t){
        return mcv;
      }    
      goto STEP_CHECK_IN_TWOLINES;
  }

STEP_CHECK_IN_TWOLINES:

  while (s != t)
  { 
    // assert(jmcv[s] > 0.1);
    // assert(jmcv[t] > 0.1);
    if(pDep[s] > pDep[t]){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }      
      s = jup[s];
    }
    else if(pDep[s] < pDep[t]){
      mcv = min(mcv, jmcv[t]);
      if(mcv == min_bound2){
        return mcv;
      }      
      t = jup[t];
    }
    else{
      mcv = min(mcv, min(jmcv[s],jmcv[t]));
      if(mcv == min_bound2){
        return mcv;
      }         
      s = jup[s];
      t = jup[t];
    }

  }

  return mcv;

}

/*
long solveMaxFlowAccVER3(long minCandi, NodePropArr np, cType s, cType t)
{
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;
  cType *paccup = np.pacc_upid;
  long *paccmcv = np.pacc_upmincv;
  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  cType depT = pDep[t];
  assert(pDep[s] >= depT);

  long mcv = MAX_LONG;
  cType prevs = 0;
  paccmcv[prevs] = MAX_LONG;

  if (pDep[s] == pDep[t])
  {
    if(paccup[s] == paccup[t]){
      goto STEP3;
    }   
    goto STEP2;
  }

  //find the first SPAN node of s and t
  cType st = t;
  if (pDep[t] % SPAN_LEN != 0)
  {
    st = paccup[t];
  }

  cType ups = paccup[s];
  depT = pDep[t];
  cType stDep = pDep[st];

  if (pDep[ups] != pDep[st])
  {
    assert(pDep[ups] > pDep[st]);

    while (pDep[ups] > stDep)
    {
      mcv = min(mcv, paccmcv[s]);
      s = ups;
      ups = paccup[ups];
      assert(pDep[ups] % SPAN_LEN == 0);
    }
  }

  assert(pDep[ups] == stDep);

  if (ups == st)
  {
    cType min_bound2 = min(paccmcv[s],paccmcv[t]);
    if(min_bound2 >= mcv || min_bound2 >=minCandi){
      return mcv;
    }

    //here: min_bound2 is the smallest currently

    assert(pDep[s] > depT);

/////////////////////////////////BRANCH 1: use min position to accelerate
    // //if know min value position, can be accelerated
    // if(pDep[s] - gpaccposmcv[s] >= depT){
    //   //if the min value of s line is achived below t
    //   if(min_bound2 == paccmcv[s]){
    //     //and the min of two lines is achieved in s line
    //     //min value must be achived before joint point of s and t, no need to search dep < depT
    //     return min_bound2;
    //   }
    //   else{
    //     //min is achieved in t line
    //     //we don't know whether min is achieved before or after joint point
    //     //DO NOTHING
    //   }
    // }
    // else{
    //   //if min of s line is achieved above t (not including t)
    //   //we don't know whether min is achieved before or after joint point

    // }
////////////////////////////////////////END of branch 1

//TODO: if we know the nearest joint dep, we can accelerate

    while (pDep[s] > depT)
    {
      mcv = min(mcv, pCV[s]);
      if(mcv == min_bound2){
        return mcv;
      }
      s = pFa[s];
    }

    if (s == t)
    {
      return mcv;
    }

    goto STEP3;
  }
  else
  {
    mcv = min(mcv, paccmcv[s]);
    if (st != t)
    {
      mcv = min(mcv, paccmcv[t]);
    }
    t = st;
    s = paccup[s];
    if(paccup[s] == paccup[t]){
      goto STEP3;
    }
    goto STEP2;
  }

  // if (depT % SPAN_LEN == 0) // t is SPAN node
  // {
  //   //it does not matter s is or not SPAN node
  //   if (pDep[s] == depT)
  //   {
  //   }
  //   else
  //   {
  //     //pDep[s] > depT
  //     prevs = 0;
  //     while (pDep[s] > depT)
  //     {
  //       mcvassert(pFa, pCV, paccup, paccmcv, s);
  //       mcv = min(mcv, paccmcv[prevs]);
  //       prevs = s;
  //       s = paccup[s];
  //       assert(pDep[s] % SPAN_LEN == 0);
  //     }

  //     assert(pDep[s] == depT);
  //     mcv = min(mcv, paccmcv[prevs]);
  //     if (s == t)
  //     {
  //       return mcv;
  //     }
  //   }
  //   //catch s
  //   assert((pDep[s] == pDep[t]) && s != t && pDep[s] % SPAN_LEN == 0);
  // }
  // else
  // {
  //   assert(s != t);
  //   if (pDep[s] == depT)
  //   {
  //     //both are not span nodes
  //     if (paccup[s] == paccup[t])
  //     {
  //       //search and return
  //       long min_bound = min(paccmcv[s], paccmcv[t]);
  //       while (s != t)
  //       {
  //         mcv = min(mcv, min(pCV[s], pCV[t]));
  //         if (mcv == min_bound)
  //         {
  //           return mcv;
  //         }
  //         s = pFa[s];
  //         t = pFa[t];
  //       }
  //       return mcv;
  //     }
  //     else
  //     {
  //       mcv = min(mcv, min(paccmcv[s], paccmcv[t]));
  //       assert((pDep[s] == pDep[t]) && s != t && pDep[s] % SPAN_LEN == 0);
  //     }
  //   }
  //   else
  //   { //pDep[s] > depT, T is not SPAN node, s may or may not SPAN node
  //     prevs = 0;
  //     while (pDep[s] > depT)
  //     {
  //       mcvassert(pFa, pCV, paccup, paccmcv, s);
  //       mcv = min(mcv, paccmcv[prevs]);
  //       prevs = s;
  //       s = paccup[s];
  //       assert(pDep[s] % SPAN_LEN == 0);
  //     }

  //     assert(pDep[s] < depT);

  //     mcv = min(mcv, min(paccmcv[prevs], paccmcv[t]));
  //     t = paccup[t];
  //     if (t == s)
  //     {
  //       return mcv;
  //     }
  //     else
  //     {
  //       assert((pDep[s] == pDep[t]) && s != t && pDep[s] % SPAN_LEN == 0); //ok
  //     }
  //   }
  // }

  // assert((pDep[s] == pDep[t]) && s != t && pDep[s] % SPAN_LEN == 0);
  // depT = MAX_LONG;

//if s and t in same depth, no matter whether SPAN node
STEP2:
  //(2) approach before the common SPAN ode
  //it means the common node is not before the next SPAN_LEN node
  assert((s != t) && (pDep[s] == pDep[t]));
  assert(paccup[s]!=paccup[t]);
  prevs = 0;
  cType prevt = 0;
  paccmcv[0] = MAX_LONG;
  while (s != t)
  {
    mcv = min(mcv, min(paccmcv[prevs], paccmcv[prevt]));
    mcvassert(pFa, pCV, paccup, paccmcv, s);
    mcvassert(pFa, pCV, paccup, paccmcv, t);

    prevs = s;
    s = paccup[s];

    prevt = t;
    t = paccup[t];

    assert(pDep[s] % SPAN_LEN == 0);
    assert(pDep[t] == pDep[s]);
  }

  //here s==t,but paccmcv[prevs] and prevt is not included in mcv

  //here

  //(3) apprach the common SPAN_LEN node, and s==t
  s = prevs;
  t = prevt;

STEP3:
  assert((s != t) && (pDep[s] == pDep[t]));
  assert(paccup[s] == paccup[t]);
  //see whether need to search
  long min_bound = min(paccmcv[s], paccmcv[t]);
  if(min_bound >= mcv || min_bound >=minCandi){
    return mcv;
  }

  //todo: if know min value position, can be accelerated

  while (s != t)
  {
    assert(pDep[s] == pDep[t]);
    mcv = min(mcv, min(pCV[s], pCV[t]));
    if (mcv == min_bound)
    {
      return mcv;
    }
    s = pFa[s];
    t = pFa[t];
  }

  return mcv;
}
*/


// long solveMaxFlowAccVER2(NodePropArr np, cType s, cType t)
// {
//   cType * pDep = np.pdep;
//   long* pCV = np.pcv;
//   cType * pFa = np.pfa;
//   cType *paccup = np.pacc_upid;
//   long *paccmcv = np.pacc_upmincv;

//   if(pDep[s] < pDep[t]){
//     cType tmp = s;
//     s = t;
//     t = tmp;
//   }

//   assert(pDep[s] >= pDep[t]);

//   long mcv = MAX_LONG;
//   cType depT = pDep[t];

//   //(1)approach the SPAN_LEN just below t
//   cType prevs = 0;
//   paccmcv[0] = MAX_LONG;
//   while(pDep[s] > depT){
//     mcv = min(mcv,paccmcv[prevs]);
//     prevs = s;
//     s = paccup[s];
//     assert(pDep[s] % SPAN_LEN == 0);
//   }

//   if(pDep[s] == depT){
//     mcv = min(mcv,paccmcv[prevs]);
//     if(s == t){
//       return mcv;
//     }
//     else{
//       //s an t same depth, and both SPAN_LEN nodes, mcv is updated
//       assert((pDep[s] == pDep[t]) && s!=t && pDep[s] %SPAN_LEN == 0);
//     }
//   }
//   else{
//     mcv = min(mcv, paccmcv[t]);
//     t = paccup[t];
//     if( t == s){
//       return mcv;
//     }
//     else{
//       //s and t same depth,and both SPAN_LEN nodes
//       assert((pDep[s] == pDep[t]) && s!=t && pDep[s] %SPAN_LEN == 0);
//     }
//   }

//   //(2) approach before the common SPAN_LEN node
//     //it means the common node is not before the next SPAN_LEN node
//   prevs = 0;
//   cType prevt = 0;
//   while(s!= t){
//     mcv = min(mcv, min(paccmcv[prevs],paccmcv[prevt]));

//     prevs = s;
//     s = paccup[s];

//     prevt = t;
//     t = paccup[t];

//     assert(pDep[s] % SPAN_LEN == 0);
//     assert(pDep[t] == pDep[s]);
//   }

//   //here

//   //(3) apprach the common SPAN_LEN node, and s==t
//   s = prevs;
//   t = prevt;
//   long min_bound = min(paccmcv[s],paccmcv[t]);
//   while(s!=t){
//     assert(pDep[s] == pDep[t]);
//     mcv = min(mcv, min(pCV[s],pCV[t]));
//     if(mcv == min_bound){
//       return mcv;
//     }
//     s = pFa[s];
//     t = pFa[t];
//   }

//   return mcv;

// }

long solveMaxFlowAccVER1(NodePropArr np, cType s, cType t)
{
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;
  cType *paccup = np.pacc_upid;
  long *paccmcv = np.pacc_upmincv;

  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  assert(pDep[s] >= pDep[t]);

  long mcv = MAX_LONG;
  cType depT = pDep[t];

  //(1)approach the nearest SPAN_LEN depth node before Depth of t (not t)
  cType prevs = s;
  while (pDep[s] > depT)
  {
    if (prevs != s)
    {
      mcv = min(mcv, paccmcv[prevs]);
    }

    prevs = s;
    s = paccup[s];
    assert(pDep[s] % SPAN_LEN == 0);
  }
  s = prevs;

  //(2) apprach common depth with T

  cType dist = pDep[s] - depT;
  // long min_bound = min(paccmcv[s],paccmcv[t]);
  while (dist-- > 0)
  {
    mcv = min(mcv, pCV[s]);
    s = pFa[s];
  }

  //(3) apprach the neartest SPAN_LEN node before common ancestor of s and t
  //here dep s == dep t but cv[s] is not included in mcv
  assert(pDep[s] == pDep[t]);
  prevs = s;
  cType prevt = t;
  while (s != t)
  {
    if (prevs != s)
    {
      assert(prevt != t);
      mcv = min(mcv, min(paccmcv[prevs], paccmcv[prevt]));
    }

    prevs = s;
    prevt = t;
    s = paccup[s];
    t = paccup[t];

    assert(pDep[s] % SPAN_LEN == 0);
    assert(pDep[t] % SPAN_LEN == 0);
  }

  //(4)approach the common ancestor of s and t
  s = prevs;
  t = prevt;

  while (s != t)
  {
    assert(pDep[s] == pDep[t]);
    mcv = min(mcv, min(pCV[s], pCV[t]));
    s = pFa[s];
    t = pFa[t];
  }

  return mcv;
}

long solveMaxFlow(NodePropArr np, cType s, cType t)
{
  long minfv = MAX_LONG;
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;

  cType depT = pDep[t];
  while (pDep[s] > depT)
  {
    // printf("pdep s %ld > dept %ld\n ",pDep[s], depT);
    // if(minfv > pCV[s]){
    //   minfv = pCV[s];
    // }
    minfv = min(minfv, pCV[s]);

    // if(minfv == 0){
    //   printf("now minv of s %ld is 0\n",s);
    //   exit(1);
    // }
    assert(pDep[s] >= 0 && pDep[s] <= N);
    assert(pDep[pFa[s]] == pDep[s] - 1);
    s = pFa[s];
  }

  cType depS = pDep[s];
  while (depS < pDep[t])
  {
    // if(minfv > pCV[t]){
    //   minfv = pCV[t];
    // }

    minfv = min(minfv, pCV[t]);

    // if(minfv == 0){
    //   printf("now minv of t %ld is 0\n",t);
    //   exit(1);
    // }
    assert(pDep[t] >= 0 && pDep[t] <= N);
    assert(pDep[pFa[t]] == pDep[t] - 1);
    t = pFa[t];
  }

  while (s != t)
  {
    // printf("s %d t %d\n",s,t);
    // long tmpv = min(pCV[s],pCV[t]);
    // if(minfv > tmpv) minfv = tmpv;

    minfv = min(minfv, min(pCV[s], pCV[t]));
    // if(minfv == 0){
    //   printf("now minv of s %ld  t %ld is 0,cv is %ld %ld, dep is %ld %ld\n",s,t,pCV[s],pCV[t],pDep[s],pDep[t]);
    //   for(int i=1; i<=8;i++){
    //     printf("n %ld dep %ld \n",i,pDep[i]);
    //   }
    //   exit(1);
    // }
    assert(pDep[s] >= 0 && pDep[s] <= N);
    assert(pDep[t] >= 0 && pDep[t] <= N);
    assert(pDep[pFa[s]] == pDep[s] - 1);
    assert(pDep[pFa[t]] == pDep[t] - 1);

    s = pFa[s];
    t = pFa[t];
  }

  return minfv;
}
// PassRes* extractResAndReset(nodeP* nodes, int N){

//   cType * nlist = walloc(N,sizeof(cType));
//   long* clist = walloc(N,sizeof(long));
//   nodeP* np;

//   for(int i=1; i<=N;i++){
//     np = nodes+i;
//     nlist[i] = np->fa;
//     clist[i] = np->cv;
//     np->s = np->cv = np->dep = np->fa = 0;
//   }

//   PassRes* res = walloc ( 1, sizeof(PassRes) );
//   res->faList = nlist;
//   res->cvList = clist;
//   return res;

// }

/*
the input condtion: one edge only appear once: 1 2 /2 1 only once as 1 2 or 2 1
node 0 not used

*/

//gen a arr with shuffled 0 - num-1
int *genShffule(int num)
{
  int *arr = walloc(num, sizeof(int));
  for (int i = 0; i < num; i++)
  {
    arr[i] = i;
  }

  for (int i = 0; i < num; i++)
  {
    int idx = i + rand() % (num - i);
    if (idx != i)
    {
      int tmpt = arr[i];
      arr[i] = arr[idx];
      arr[idx] = tmpt;
    }
  }

  return arr;
}

int main(argc, argv)

    int argc;
char *argv[];

{
  //PARAM
  int P = 30; // P% percent total passes in mode 1, the remaining in mode 2
  int total = 500; //number of total passes
  long numOfPairs = 100; //number of random (s,t) pairs
  //PARAM END

  if (argc > 2)
  {
    printf("Usage: %s [update frequency]\n", argv[0]);
    exit(1);
  }

  printf("c\nc hi_treem version 0.9\n");
  printf("c Copyright C by nsyncw, nsyncw@126.com\nc\n");

  parse(&N, &M, &nodes);

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", N, M);
  initrand();
  //init SPAN_LEN
  SPAN_LEN = (int)(sqrt(N));
  // nodeP* np = nodes;
  //init results:

  // nodeP* np = nodes;
  // for(int i=1; i<=N;i++){
  //   np = nodes+i;
  //   if(np->edges != NULL){
  //     // printf("order for %d before\n",i);
  //     deOrderEdgeByRandomCap(np);
  //     // printf("---------%ld: %ld %ld %ld %ld\n",i,np->s,np->cv,np->dep,np->fa);
  //     // printf("order for %d after\n",i);
  //     // orderEdgeByCap(np->edges);
  //     for(int j=0; j<np->nIdx; j++){
  //         // if(eh->endNode < i){
  //         //   continue;
  //         // }
  //         edgeP* eh = np->edges+ np->orderedEdges[j];
  //         printf("a %ld %ld %ld\n", i,eh->endNode, eh->cap);
  //         // printf("   --------   %ld %ld %ld \n",eh->endNode,eh->rev->endNode, eh->rev->cap);
  //     }
  //   }
  // }

  // exit(0);

  assert(total <= MAX_PASSES);

  cType *roots = (cType *)walloc(total + 2, sizeof(cType));
  assert(roots != NULL);
  for (int i = 0; i < total; i++)
  {
    allResults[i].pfa = (cType *)walloc(N + 2, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(N + 2, sizeof(cType));
    allResults[i].pcv = (long *)walloc(N + 2, sizeof(long));
    allResults[i].ps = (short *)walloc(N + 2, sizeof(short));
    allResults[i].pacc_upid = (cType *)walloc(N + 2, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)walloc(N + 2, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)walloc(N + 2, sizeof(cType));
    allResults[i].pacc_jointid = (cType *)walloc(N + 2, sizeof(cType));
    allResults[i].pacc_jointmcv = (long *)walloc(N + 2, sizeof(long));


    // assert(allResults[i].pfa != NULL);
    // assert(allResults[i].pdep != NULL);
    // assert(allResults[i].pcv != NULL);
    // assert(allResults[i].ps != NULL);

    memset(allResults[i].pfa, 0, (N + 2) * sizeof(cType));
    memset(allResults[i].pdep, 0, (N + 2) * sizeof(cType));
    memset(allResults[i].pcv, 0, (N + 2) * sizeof(long));
    memset(allResults[i].ps, 0, (N + 2) * sizeof(short));
    memset(allResults[i].pacc_upid, 0, (N + 2) * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, (N + 2) * sizeof(long));
    memset(allResults[i].pacc_pos_upmincv, 0, (N + 2) * sizeof(cType));
    memset(allResults[i].pacc_jointid, 0, (N + 2) * sizeof(cType));
    memset(allResults[i].pacc_jointmcv, 0, (N + 2) * sizeof(long));
  }

  // printf("debug before\n");
  // for(int i=0; i<100;i++){
  //   cType root = 1+(mrand()%N);
  //   printf("idx %ld rand %lu\n",randNumIdx, root);
  // }

  // exit(0);

  double tm;
  double totalProcTime = 0;
  //calculate total cap of one node
  preMarkCut(1 + ((mrand() * mrand()) % N));
  cType root = 1 + ((mrand() * mrand()) % N);
  printf("c root for all passes is %ld\n",root);
  for (int ipass = 0; ipass < total; ipass++)
  {
    // printf("the %d times\n",i);
    gpfa = allResults[ipass].pfa;
    gpdep = allResults[ipass].pdep;
    gpcv = allResults[ipass].pcv;
    gps = allResults[ipass].ps;
    gpaccup = allResults[ipass].pacc_upid;
    gpaccmcv = allResults[ipass].pacc_upmincv;
    gpaccposmcv = allResults[ipass].pacc_pos_upmincv;
    gpaccjointid = allResults[ipass].pacc_jointid;
    gpaccjointmcv = allResults[ipass].pacc_jointmcv;

    mode = ipass < P * total / 100 ? 1 : 2;
    initrand();
    
    roots[ipass] = root;
    gpdep[root] = 0;
    // printf("pass %d, randidx %ld, root is %ld\n",i, randNumIdx,root);
    // printf("root fa %ld\n",allResults[i].pfa[root]);
    tm = timer();
    markCut(root);
    gpcv[root] = MAX_LONG;
    // gpaccjointid[root] = 2;//just mark root as a joint node
    buildAcc(root, root, MAX_LONG,gpdep[root],root,MAX_LONG);

    // printf("---\n");
    // for(int k=0; k<=N; k++){
    //   printf("pass0 pFa %d is %ld\n",k,allResults[0].pfa[k]);
    // }

    totalProcTime += timer() - tm;
    printf("proctime for onepass: %10.06f\n", timer() - tm);
    // printf("a root fa %ld\n",allResults[i].pfa[root]);
    if (ipass % 10 == 0)
    {
      printf("the %d passes\n", ipass);
    }
  }

  printf("_TAR_X: preprocess times %10.6f\n", totalProcTime);

  gpfa = NULL;
  gpdep = NULL;
  gpcv = NULL;
  gps = NULL;
  gpaccup = NULL;
  gpaccmcv = NULL;

  //------FOR DEBUG
  for (int i = 0; i < total; i++)
  {
    // printf("------------------for pass %d, root is %ld\n",i,roots[i]);
    if (allResults[i].pfa[roots[i]] != 0)
    {
      printf("root %ld 's fa is not 0,is %ld", roots[i], allResults[i].pfa[roots[i]]);
      // for(int k=0; k<=N; k++){
      //   printf("pFa %d is %ld\n",k,allResults[i].pfa[k]);
      // }
      exit(1);
    }
    for (int j = 1; j <= N; j++)
    {
      assert(!(allResults[i].pcv[j] == 0 && allResults[i].pdep[j] != 0));
      // if(allResults[i].pcv[j] == 0 && allResults[i].pdep[j] != 0){
      //   printf("root is %ld node %ld cv is 0, s %ld fa %ld dep %ld\n",roots[i],j,allResults[i].ps[j],allResults[i].pfa[j],allResults[i].pdep[j]);
      //   exit(1);
      // }

      // printf("node %ld is %ld fa%ld %ld\n",j,allResults[i].ps[j],allResults[i].pfa[j],allResults[i].pdep[j]);
    }
  }
  //-------FOR DEBUG END

  // exit(0);
  initrand();
  // for(int i=0; i<100; i++){
  //   printf("idx %ld rand %lu \n",randNumIdx, mrand());
  // }

  double totalTime = 0;
  long mv = MAX_LONG;

  double curTime = 0;
  cType ns, nt;

  // int* arr_total_index = genShffule(total);

  for (int ipair = 0; ipair < numOfPairs;)
  {
    // printf("%d\n",i);
    ns = 1 + ((mrand() * mrand()) % (N));
    nt = 1 + ((mrand() * mrand()) % (N));
    if (ns != nt)
    {
      mv = MAX_LONG;
      curTime = timer();
      for (int j = 0; j < total; j++)
      {
        long tmp = solveMaxFlowAccVER4(mv, allResults[j], ns, nt);
        if (mv > tmp)
        {
          mv = tmp;
        }
      }
      curTime = timer() - curTime;
      totalTime += curTime;
      ipair++;
      printf("hi_treem_res(n,s,mflow,tm) %lu %lu %12.04f %12.06f1\n", ns, nt, 1.0 * mv, curTime);
      // printf("echo \"c flow:\t%lu\t%12.04f --treem pair %d/%d \"\n",mv,curTime,ipair,numOfPairs);
      // printf("bash solven100000.sh %lu %lu | grep -E \"flow:\"\n",ns,nt);
      // printf("time is %10.6f, maxflow is %ld for %ld %ld\n",curTime,mv,ns,nt);
    }
  }

  printf("#run ok! average time %10.6f\n", totalTime / numOfPairs);
  // cc = allocDS();
  // if ( cc ) { fprintf ( stderr, "Allocation error\n"); exit ( 1 ); }

  // init();

  // forAllNodes(i) {
  //   printf("for node %ld: \n",i - nodes + 1);

  //   forAllArcs(i,a){
  //     if(a->rev->resCap > 0){
  //       printf("\t\t%ld %ld %ld \n",a->rev->head - nodes +1 , a->head - nodes +1,a->rev->resCap);
  //     }
  //   }
  // }

  exit(0);

  //   printf("c init done\n");
  //   stageOne ( );

  //   t2 = timer() - t2;

  //   printf ("c flow:       %12.01f\n", flow);

  // #ifndef CUT_ONLY
  //   stageTwo ( );

  //   t = timer() - t;

  //   printf ("c time:        %10.2f\n", t);

  // #endif

  //   printf ("c cut tm:      %10.2f\n", t2);

  // #ifdef CHECK_SOLUTION

  //   /* check if you have a flow (pseudoflow) */
  //   /* check arc flows */
  //   forAllNodes(i) {
  //     forAllArcs(i,a) {
  //       if (cap[a - arcs] > 0) /* original arc */
  // 	if ((a->resCap + a->rev->resCap != cap[a - arcs])
  // 	    || (a->resCap < 0)
  // 	    || (a->rev->resCap < 0)) {
  // 	  printf("ERROR: bad arc flow\n");
  // 	  exit(2);
  // 	}
  //     }
  //   }

  //   /* check conservation */
  //   forAllNodes(i)
  //     if ((i != source) && (i != sink)) {
  // #ifdef CUT_ONLY
  //       if (i->excess < 0) {
  // 	printf("ERROR: nonzero node excess\n");
  // 	exit(2);
  //       }
  // #else
  //       if (i->excess != 0) {
  // 	printf("ERROR: nonzero node excess\n");
  // 	exit(2);
  //       }
  // #endif

  //       sum = 0;
  //       forAllArcs(i,a) {
  // 	if (cap[a - arcs] > 0) /* original arc */
  // 	  sum -= cap[a - arcs] - a->resCap;
  // 	else
  // 	  sum += a->resCap;
  //       }

  //       if (i->excess != sum) {
  // 	printf("ERROR: conservation constraint violated\n");
  // 	exit(2);
  //       }
  //     }

  //   /* check if mincut is saturated */
  //   aMax = dMax = 0;
  //   for (l = buckets; l < buckets + n; l++) {
  //     l->firstActive = sentinelNode;
  //     l->firstInactive = sentinelNode;
  //   }
  //   globalUpdate();
  //   if (source->d < n) {
  //     printf("ERROR: the solution is not optimal\n");
  //     exit(2);
  //   }

  //   printf("c\nc Solution checks (feasible and optimal)\nc\n");
  // #endif

  // #ifdef PRINT_STAT
  //     printf ("c pushes:      %10ld\n", pushCnt);
  //     printf ("c relabels:    %10ld\n", relabelCnt);
  //     printf ("c updates:     %10ld\n", updateCnt);
  //     printf ("c gaps:        %10ld\n", gapCnt);
  //     printf ("c gap nodes:   %10ld\n", gNodeCnt);
  //     printf ("c\n");
  // #endif

  // #ifdef PRINT_FLOW
  //     printf ("c flow values\n");
  //     forAllNodes(i) {
  //       ni = nNode(i);
  //       forAllArcs(i,a) {
  // 	na = nArc(a);
  // 	if ( cap[na] > 0 )
  // 	  printf ( "f %7ld %7ld %12ld\n",
  // 		  ni, nNode( a -> head ), cap[na] - ( a -> resCap )
  // 		  );
  //       }
  //     }
  //     printf("c\n");
  // #endif

  // #ifdef PRINT_CUT
  //   globalUpdate();
  //   printf ("c nodes on the sink side\n");
  //   forAllNodes(j)
  //     if (j->d < n)
  //       printf("c %ld\n", nNode(j));

  // #endif

  // exit(0);
}
