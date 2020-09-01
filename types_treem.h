/* defs.h */

//#ifdef EXCESS_TYPE_LONG
///typedef unsigned long excessType;
//#else
typedef unsigned long long int excessType; /* change to double if not supported */
//#endif

typedef unsigned long cType;
typedef unsigned int sType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap;
   cType w; /*weight*/
   cType avgCV; /*average CV of cut set where this edge belongs to */
   long tmp;

   struct edgeProp* rev; //reverse

}edgeP;


typedef  /* node */
   struct nodeProp
{
   edgeP* edges;
   cType maxEdges;
   cType nIdx;
   cType totalCap;
   
   sType* orderedEdges;

} nodeP;

typedef
struct NodePropExtra_
{
   cType fa;
   cType dep;
   long cv;
   short s;


} NodePropExtra;


typedef
struct NodePropArr_{

   cType * pfa;
   cType * pdep;
   long * pcv;
   short * ps;
   cType * pacc_upid; // up node id
   long* pacc_upmincv; // min of [currnet node, upnode)
   cType * pacc_pos_upmincv; // the diff of depth of nearest min value
   cType * pacc_jointid; //the nearest joint id;
   long* pacc_jointmcv; //the min of [cur, joint node)
   
} NodePropArr;


// typedef  /* arc */
//    struct arcSt
// {
//    cType           resCap;          /* residual capasity */
//    struct nodeSt   *head;           /* arc head */
//    struct arcSt    *rev;            /* reverse arc */
// }
//   arc;




// typedef /* bucket */
//    struct bucketSt
// {
//   node             *firstActive;      /* first node with positive excess */
//   node             *firstInactive;    /* first node with zero excess */
// } bucket;

