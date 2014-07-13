#include<stdio.h>
#include<math.h>
#include<cmath>
#include<string.h>
//#include "ran2.c"
#define pi 3.141592654
#define noli 10                              /* num of lines on one face of bipyramid */                                
#define blen 1.3000                          /* bond length to make bipyramid         */
//#define nver 3*(noli*noli)+2                 /* num of vertices on the surface        */
//#define ntr (nver-2)*2                       /* number of triangles                   */
#define ttether (nver-2)*3                     /* number of tethers                       */
#define conf_int 1000                         /* to print configuration data file      */
#define itime 1                              /* initial point of simulation           */
#define ftime 10000                           /* final point of simulation             */
#define kappa 10.0                           /* bending rigidity coefficient          */
#define pr 0.0                               /* pressure difference between in & out  */
#define fangle -0.5000                       /* minimum face angle allowed            */
#define beta 1.0000                                                               
                                             /* data structure for vertices, triangles & tethers */
using namespace std;

int nver;
int ntr;

typedef struct vertex { int nonei,vneipt[10],vneitr[10];
                        double vcoord[4],mcur,cur1,cur2,totarea,vnor[4];
                      } vertex;
typedef struct triangle { int vert[4],li[4];
                          double fnor[4],ar,vol;
                        }triangle;
typedef struct tether { int tr,sep[3];
                   } tether;
typedef struct tetherm { int tr,sep[3];
                    } tetherm;  

void STARTUP(vertex *,triangle *,tether *,tetherm * );                      /* initialisation with bipyramid     */
void JAVAVIEWDAT(int , vertex *,triangle *,tether *,tetherm *);             /* configuration data fiel           */
void CURVCALC(int ,vertex *,triangle *);                                /* curvature calculation at vertices */
void MOVEVERTEX(int ,vertex *,triangle *,tether *,tetherm *);               /* to move a vertex                  */
void FLIPPING(int ,vertex *,triangle *,tether *,tetherm *);                 /* to flip a tether                    */
void areacal(int ,triangle *,vertex *);                                 /* to calculate area of a triangle   */
void onlyarea(int ,triangle *,vertex *);
int faceangchk(int ,triangle *,tether *,tetherm *);                         /* to check minm face angle          */
void system_energy(int ,triangle *,vertex *);
float ran2(long *idum);

FILE *fp1;

long int seed= -637648728;

/*                                           MAIN PROGRAM STARTS HERE                                         */
/*˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚*/
void mainfunc()
  {
        int i,mcs,rand,iloop;
        float randf;
        fp1=fopen("mem_prop.dat","w");
        vertex ver[(nver+1)*3]; vertex (*vptr)[(nver+1)*3];      
        triangle tri[ntr+1];    triangle (*trptr)[ntr+1];
        tether lin[ttether+1];      tether (*liptr)[ttether+1];
        tetherm linm[ttether+1];    tetherm (*limptr)[ttether+1];
        vptr=&ver; trptr=&tri; liptr=&lin; limptr=&linm; 
        
        mcs=0;
        STARTUP(ver,tri,lin,linm);
        for(i=1;i<=nver;i++) CURVCALC(i,ver,tri);
        JAVAVIEWDAT(mcs/conf_int,ver,tri,lin,linm);


        mcsloop:
        for(mcs=itime;mcs<=ftime;mcs++){
        innerloop:
        for(iloop=1;iloop<=ttether;iloop++){
             rand=rint(ran2(&seed)*nver);
             if(rand>0 && rand <=nver) MOVEVERTEX(rand,ver,tri,lin,linm); 
             rand=rint((1-2*ran2(&seed))*ttether);
             if((rand!=0) && (abs(rand) <=ttether)) FLIPPING(rand,ver,tri,lin,linm);
            }
         if(mcs%conf_int==0){ JAVAVIEWDAT(mcs/conf_int,ver,tri,lin,linm);
                              system_energy(mcs,tri,ver); }
         }
  }

/*                                  MAIN PROGRAM ENDS HERE                                                       */
/*˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚˚*/  

/*---------------------------------------------------------------------------------------------------------------*/              
/*                          FUNCTION TO MAKE A TRIANGULATED BIPYRAMID                                            */
/*---------------------------------------------------------------------------------------------------------------*/

void STARTUP(vertex *vptr,triangle *trptr,tether *liptr,tetherm *limptr )
  {
       int N,lz,ll,ls,ipn,nt,nt2,i,j,np,in,ite,nbm[nver+1],ip[nver+1][10],tl[nver+1][nver+1];
       int trl[nver+1][nver+1],nlin,k,n1,n2,trno;
       double dl,H,theta,b;
       double x[nver+1],y[nver+1],z[nver+1];
       N=noli;b=blen;
       ipn=0;
       nt=3*N*(N+1)/2+1;
       dl=N*b;
       theta=asin(1.0/sqrt(3.0));
       H=dl*cos(theta);
                                                   
       up_half_face1:for(lz=0;lz<=N;lz++){                                    /* points on first face of upper   */
                      if(lz<1) ll=0;                                          /* half plane                      */
                      if(lz>=1)ll=lz-1;
                      for(ls=0;ls<=ll;ls++){
                             ipn=ipn+1;
                             z[ipn]=H-b*lz*cos(theta);
                             x[ipn]=b*lz*sin(theta)/2.0;
                             y[ipn]=-dl/2.0+ls*b+(N-lz)*b/2.0;}}
      
       up_half_face2:for(lz=1;lz<=N;lz++){                                    /*points on second face of upper    */
                      for(ls=0;ls<=lz-1;ls++){                                /* half plane                       */
                             ipn=ipn+1;
                             z[ipn]=H-b*lz*cos(theta);
                             x[ipn]=b*(lz-ls)*sqrt(3.0)/2.0-b*lz*sin(theta);
                             y[ipn]=b*lz/2.0-ls*b/2.0; }}
 
       up_half_face3:for(lz=1;lz<=N;lz++){                                    /* points on third face of upper    */
                     for(ls=0;ls<=lz-1;ls++){                                 /* half plane                       */
                            ipn=ipn+1;
                            z[ipn]=H-b*lz*cos(theta);
                            x[ipn]=b*ls*sqrt(3.0)/2.0-b*lz*sin(theta);
                            y[ipn]=-ls*b/2.0;}}

       nbm[1]=3;                                                              /* number of neighbors for top &    */
       nbm[nt+1]=3;                                                           /* bottom particle is 3             */

       for(i=1;i<=nver;i++){
       for(j=1;j<=9;j++)ip[i][j]=0;}

       up_nei_loop1:for(i=1;i<=3;i++){                                       /* neighbhors of particle 1 is fixed */
                   ip[1][i]=(i-1)*N*(N+1)/2+2;                               /* in every loop (2,17,32 for l=5 )  */
                   lz=1;ls=0;
                   np=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+2;                      /* the point nearest to the top point*/
                   nbm[np]=6;                                                /* all of these have 6 neigh each    */
                   in=i+1;  
                   if(in>3) in=1;
                   ip[np][1]=(in-1)*N*(N+1)/2+2;
                   ip[np][2]=1;
                   in=i-1;
                   if(in<1)in=3;

                   ip[np][3]=(in-1)*N*(N+1)/2+2;
                   ip[np][4]=(in-1)*N*(N+1)/2+2+2;
                   ip[np][5]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+2;
                   ip[np][6]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2;
          
            up_nei_loop2:for(lz=2;lz<=N;lz++){
            up_nei_loop3:for(ls=0;ls<=lz-1;ls++){
                         np=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+2;
                         nbm[np]=6;

                         if(ls!=0 && lz!=N){
                              if(ls!=lz-1)
                                   {ip[np][1]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2;
                                   ip[np][2]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2;}
                              else
                                   {ite=i*N*(N+1)/2+lz*(lz-1)/2+2;
                                   if(ite > nt)ite=ite-nt+1;
                                   ip[np][1]=ite;
                                   ite=i*N*(N+1)/2+(lz-1)*(lz-2)/2+2;
                                   if(ite > nt)ite=ite-nt+1;
                                   ip[np][2]=ite; }

                              ip[np][3]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls-1+2;
                              ip[np][4]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls-1+2;
                              ip[np][5]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+2;
                              ip[np][6]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2;}
                         
                        else if(ls==0 && lz!=N){
                              in=i-1;
                              if(in < 1)in=3;

                              ip[np][1]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2;
                              ip[np][2]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2;
                              ip[np][3]=(in-1)*N*(N+1)/2+lz*(lz-1)/2+lz-1+2;
                              ip[np][4]=(in-1)*N*(N+1)/2+lz*(lz+1)/2+lz+2;
                              ip[np][5]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+2;
                              ip[np][6]=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2; }

                       else if(ls!=0 && lz==N){
                              if(ls != lz-1){
                                   ip[np][1]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2;
                                   ip[np][2]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2;
                                   ip[np][6]=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+ls+2+nt;}
                              else {ite=i*N*(N+1)/2+lz*(lz-1)/2+2;
                                   if(ite > nt)ite=ite-nt+1;
                                   ip[np][1]=ite;
                                   ite=i*N*(N+1)/2+(lz-1)*(lz-2)/2+2;
                                   if(ite > nt)ite=ite-nt+1;
                                   ip[np][2]=ite;
                                   in=i+1;
                                   if(in>3)in=1;
                                   ip[np][6]=(in-1)*N*(N-1)/2+(N-1)*(N-2)/2+2+nt; }

                              ip[np][3]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls-1+2;
                              ip[np][4]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls-1+2;
                              ip[np][5]=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+ls-1+2+nt;}

                     else if(ls==0 && lz==N)
                            { in=i-1;
                              if(in<1)in=3;
                              ip[np][1]=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2;
                              ip[np][2]=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2;
                              ip[np][3]=(in-1)*N*(N+1)/2+lz*(lz-1)/2+lz-1+2;
                              ip[np][4]=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+2+nt;
                              nbm[np]=4; }          
          }}} 

          low_half_face1:for(lz=0;lz<=N-1;lz++){                               /*  points in the lower half faces       */
                         if(lz < 1)ll=0;
                         if(lz > 1)ll=lz-1;
                         for(ls=0;ls<=ll;ls++){
                                 ipn=ipn+1;
                                 z[ipn]=-H+b*lz*cos(theta);
                                 x[ipn]=b*lz*sin(theta)/2.0;
                                 y[ipn]=-dl/2.0+ls*b+(N-lz)*b/2.0;}}
          low_half_face2:for(lz=1;lz<=N-1;lz++){
                         for (ls=0;ls<=lz-1;ls++){
                                 ipn=ipn+1;
                                 z[ipn]=-H+b*lz*cos(theta);
                                 x[ipn]=b*(lz-ls)*sqrt(3.0)/2.0-b*lz*sin(theta);
                                 y[ipn]=b*lz/2.0-ls*b/2.0; }}
          low_half_face3:for (lz=1;lz<=N-1;lz++){
                         for(ls=0;ls<=lz-1;ls++){
                                 ipn=ipn+1;
                                 z[ipn]=-H+b*lz*cos(theta);
                                 x[ipn]=b*ls*sqrt(3.0)/2.0-b*lz*sin(theta);
                                 y[ipn]=-ls*b/2.0;}}

          lnei_loop1: for(i=1;i<=3;i++){                                      /* neighbhors of points in lowes half faces*/
                      ip[nt+1][4-i]=(i-1)*N*(N-1)/2+2+nt;
                      lz=1 ; ls=0;
                      np=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+2+nt;
                      nbm[np]=6;

                      in=i+1;
                      if(in>3) in=1;

                      ip[np][6]=(in-1)*N*(N-1)/2+2+nt;
                      ip[np][5]=1+nt;

                      in=i-1;
                      if(in < 1) in=3;

                      ip[np][4]=(in-1)*N*(N-1)/2+2+nt;
                      ip[np][3]=(in-1)*N*(N-1)/2+4+nt;
                      ip[np][2]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+2+nt;
                      ip[np][1]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt;
                      nt2=3*N*(N-1)/2+1;

                lnei_loop2:for(lz=2;lz<=N-1;lz++){
                lnei_loop3:for(ls=0;ls<=lz-1;ls++){
                           np=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+2+nt;
                           nbm[np]=6;

                           if(ls!=0 && lz!= N-1){
                              if(ls !=lz-1){
                                  ip[np][6]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt;
                                  ip[np][5]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt;}
                              else{
                                  ite=i*N*(N-1)/2+lz*(lz-1)/2+2;
                                  if(ite > nt2)ite=ite-nt2+1;
                                  ip[np][6]=ite+nt;
                                  ite=i*N*(N-1)/2+(lz-1)*(lz-2)/2+2;
                                  if(ite > nt2)ite=ite-nt2+1;
                                  ip[np][5]=ite+nt;}
                              ip[np][4]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+1+nt; 
                              ip[np][3]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+1+nt;
                              ip[np][2]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+2+nt;
                              ip[np][1]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt;}

                        else if(ls==0 && lz!= N-1){
                              in=i-1;
                              if(in < 1)in=3;
                              ip[np][6]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt;
                              ip[np][5]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt;
                              ip[np][4]=(in-1)*N*(N-1)/2+lz*(lz-1)/2+lz+1+nt;
                              ip[np][3]=(in-1)*N*(N-1)/2+lz*(lz+1)/2+lz+2+nt;
                              ip[np][2]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+2+nt;
                              ip[np][1]=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt;}

                        else if(ls !=0 && lz==N-1){
                              if(ls != lz-1){
                                ip[np][6]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt;
                                ip[np][5]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt;}
                              else{
                                ite=i*N*(N-1)/2+lz*(lz-1)/2+2;
                                if(ite > nt2) ite=ite-nt2+1;
                                ip[np][6]=ite+nt;
                                ite=i*N*(N-1)/2+(lz-1)*(lz-2)/2+2;
                                if(ite > nt2) ite=ite-nt2+1;
                                ip[np][5]=ite+nt; }
                             
                              ip[np][4]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+1+nt;
                              ip[np][3]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+1+nt;
                              ip[np][2]=(i-1)*N*(N+1)/2+N*(N-1)/2+ls+2;
                              ip[np][1]=(i-1)*N*(N+1)/2+N*(N-1)/2+ls+3;}

                        else if(ls==0 && lz==N-1){
                                in=i-1;
                                if(in < 1) in=3;
                                ip[np][6]=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt;
                                ip[np][5]=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt;
                                ip[np][4]=(in-1)*N*(N-1)/2+lz*(lz-1)/2+lz+1+nt;
                                ip[np][3]=(in-1)*N*(N+1)/2+N*(N-1)/2+N+1;
                                ip[np][2]=(i-1)*N*(N+1)/2+N*(N-1)/2+2;
                                ip[np][1]=(i-1)*N*(N+1)/2+N*(N-1)/2+3;}
                      }}}
     
        for(i=1;i<=nver;i++){ 
        for(j=1;j<=10;j++){(vptr+i)->vneipt[j]=0;
                            (vptr+i)->vneitr[j]=0;}}
              
        for(i=1;i<=nver;i++){                                                 /* feed vertex coord to data   struct */
               (vptr+i)->vcoord[1]=x[i];
               (vptr+i)->vcoord[2]=y[i];
               (vptr+i)->vcoord[3]=z[i];
                (vptr+i)->nonei=nbm[i];                                       /* store the num of neigh             */
                for(j=1;j<=nbm[i];j++) (vptr+i)->vneipt[j]=ip[i][j];}         /* neighboring vertices of each vertex*/
          
        for(i=1;i<=nver;i++){
        for(j=1;j<=nver;j++) { tl[i][j]=0;trl[i][j]=0;}}

        for(i=1;i<=ntr;i++){           
               (trptr+i)->vert[1]=0; 
               (trptr+i)->vert[2]=0;
               (trptr+i)->vert[3]=0;}
           
        nt=0; nlin=0;                                                        /* finding out triangles & tethers        */
        for(i=1;i<=nver;i++){
        for(j=1;j<=nbm[i];j++){                
                  k=j+1; if(j==nbm[i]) k=1;
                  n1=ip[i][j]; n2=ip[i][k];

                  if(trl[i][n1]==0)
                   { nt=nt+1;
                     (trptr+nt)->vert[1]=i; (trptr+nt)->vert[2]=n1; (trptr+nt)->vert[3]=n2;
                     trl[i][n1]=nt; trl[n1][n2]=nt; trl[n2][i]=nt;
                     (vptr+i)->vneitr[j]=nt;

                     if(tl[i][n1]==0){
                           nlin=nlin+1;
                           (liptr+nlin)->sep[1]=i;(liptr+nlin)->sep[2]=n1;
                           (limptr+nlin)->sep[1]=n1; (limptr+nlin)->sep[2]=i;
                           tl[i][n1]=nlin; tl[n1][i]=-nlin;}
                        
                     if(tl[n1][n2]==0){
                           nlin=nlin+1; 
                           (liptr+nlin)->sep[1]=n1;(liptr+nlin)->sep[2]=n2;
                           (limptr+nlin)->sep[1]=n2;(limptr+nlin)->sep[2]=n1;
                           tl[n1][n2]=nlin ; tl[n2][n1]=-nlin;}
                        
                     if(tl[n2][i]==0){
                           nlin=nlin+1;
                           (liptr+nlin)->sep[1]=n2;(liptr+nlin)->sep[2]=i;
                           (limptr+nlin)->sep[1]=i;(limptr+nlin)->sep[2]=n2;
                           tl[n2][i]=nlin ; tl[i][n2]=-nlin;} 
                     
                    if(tl[i][n1]>0) (liptr+tl[i][n1])->tr=nt ; 
                    else (limptr-tl[i][n1])->tr=nt ;
                      
                    if(tl[n1][n2]>0) (liptr+tl[n1][n2])->tr=nt ; 
                    else (limptr-tl[n1][n2])->tr=nt ; 
                    if(tl[n2][i]>0) (liptr+tl[n2][i])->tr=nt ; 
                    else (limptr-tl[n2][i])->tr=nt ;  

                    (trptr+nt)->li[1]=tl[i][n1]; 
                    (trptr+nt)->li[2]=tl[n1][n2]; 
                    (trptr+nt)->li[3]=tl[n2][i];           
                  }
                 else
                   { trno=trl[i][n1];
                     if(tl[i][n1]>0) (liptr+tl[i][n1])->tr=trno;
                     else (limptr-tl[i][n1])->tr=trno;
                     if(tl[n1][n2]>0) (liptr+tl[n1][n2])->tr=trno;
                     else (limptr-tl[n1][n2])->tr=trno;
                     if(tl[n2][i]>0) (liptr+tl[n2][i])->tr=trno;
                     else (limptr-tl[n2][i])->tr=trno;
                     (vptr+i)->vneitr[j]=trno;  }
              }}

         for(i=1;i<=nver;i++) (vptr+i)->totarea=0.00;
         for(i=1;i<=ntr;i++) (trptr+i)->ar=0.0;
         for(i=1;i<=ntr;i++)  areacal(i,trptr,vptr);                            /* find area of all triangles    */
  }
/*---------------------------------------------------------------------------------------------------------------*/
/*                          FUNCTION TO MAKE THE AREA OF A TRIANGLE                                              */
/*---------------------------------------------------------------------------------------------------------------*/

void areacal(int tr,triangle *trptr,vertex *vptr)              // passing ith triangle and its vertex is enough 
{
     int i,j,k,n;
     double r1[4],r2[4],r3[4],r21[4],r31[4],ax,ay,az,area;
     ax=0.0; ay=0.0; az=0.0;
 
     i=(trptr+tr)->vert[1];  j=(trptr+tr)->vert[2]; k=(trptr+tr)->vert[3];              /* vertices of triangle  */
     for(n=1;n<=3;n++){ r1[n]=(vptr+i)->vcoord[n];                                    /* coord of vertices     */
                  r2[n]=(vptr+j)->vcoord[n] ; 
                  r3[n]=(vptr+k)->vcoord[n]; }   
     for(n=1;n<=3;n++){ r21[n]=r2[n]-r1[n]; r31[n]=r3[n]-r1[n];}             
     (vptr+i)->totarea=(vptr+i)->totarea-(trptr+tr)->ar*0.3333333;                
     (vptr+j)->totarea=(vptr+j)->totarea-(trptr+tr)->ar*0.3333333;
     (vptr+k)->totarea=(vptr+k)->totarea-(trptr+tr)->ar*0.3333333;
     ax=(r21[2]*r31[3]-r21[3]*r31[2])*0.5;
     ay=(r21[3]*r31[1]-r21[1]*r31[3])*0.5;
     az=(r21[1]*r31[2]-r21[2]*r31[1])*0.5;
     area=sqrt(ax*ax+ay*ay+az*az);
     (trptr+tr)->ar=area; 
     ax=ax/area ; ay=ay/area ;  az=az/area;                                             /* normalised area vector */
     (trptr+tr)->fnor[1]=ax; (trptr+tr)->fnor[2]=ay; (trptr+tr)->fnor[3]=az;
     (vptr+i)->totarea=(vptr+i)->totarea+(trptr+tr)->ar*0.3333333;                
     (vptr+j)->totarea=(vptr+j)->totarea+(trptr+tr)->ar*0.3333333;
     (vptr+k)->totarea=(vptr+k)->totarea+(trptr+tr)->ar*0.3333333;
    (trptr+tr)->vol=(r1[1]*(r2[2]*r3[3]-r2[3]*r3[2])+r1[2]*(r2[3]*r3[1]
                    -r2[1]*r3[3])+r1[3]*(r2[1]*r3[2]-r2[2]*r3[1]))/6.0;
}
/*----------------------------------------------------------------------------------------------------------------*/
/*                              FUNCTION TO WRITE THE FILES FOR JAVAVIEW                                          */
/*----------------------------------------------------------------------------------------------------------------*/
void JAVAVIEWDAT(int num, vertex *vptr,triangle *trptr, tether *liptr,tetherm *limptr)
{
    int i;
    char fname[100];
    FILE *fp;

    sprintf(fname,"jv.%d.jvx",num);
    fp=fopen(fname,"w");
    fprintf(fp,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n");
    fprintf(fp,"<!DOCTYPE jvx-model SYSTEM   \"http://www.javaview.de/rsrc/jvx.dtd\">\n");
    fprintf(fp,"<jvx-model>\n");
    fprintf(fp,"<geometries>\n");
    fprintf(fp,"<geometry name=\"membrane\">\n");
    fprintf(fp,"<pointSet dim=\"3\" color=\"show\" point=\"hide\">\n");
    fprintf(fp,"<points num=\" %d \">\n",nver);

    for(i=1;i<=nver;i++)
     fprintf(fp,"<p>\t%f\t%f\t%f</p>\n",(vptr+i)->vcoord[1],
      (vptr+i)->vcoord[2],(vptr+i)->vcoord[3]);
    fprintf(fp,"<thickness>2.0</thickness>\n");
    fprintf(fp,"</points>\n");
    fprintf(fp,"<colors num=\"%d\">\n",nver);
    for(i=1;i<=nver;i++) fprintf(fp,"<c>187 247 245\t</c>\n");
    fprintf(fp,"</colors>\n");
    fprintf(fp,"</pointSet>\n");
    fprintf(fp,"<lineSet line=\"hide\">\n");
    fprintf(fp,"<lines>\n");
    for(i=1;i<=ttether;i++) fprintf(fp,"<l>\t%d\t%d</l> \n",
     (liptr+i)->sep[1]-1,(liptr+i)->sep[2]-1);
    fprintf(fp,"<thickness>2.0</thickness>\n");
    fprintf(fp,"<color> 255 0 0 </color>\n");
    fprintf(fp,"</lines>\n");
    fprintf(fp,"</lineSet>\n");
    fprintf(fp,"<faceSet face=\"show\" edge=\"hide\">\n");
    fprintf(fp,"<faces>\n");
    for(i=1;i<=ntr;i++)fprintf(fp,"<f>\t%d\t%d\t%d</f>\n",
    ((trptr+i)->vert[1])-1,((trptr+i)->vert[2])-1,((trptr+i)->vert[3])-1);
    fprintf(fp,"<color> 187  247 245 </color>\n");
    fprintf(fp,"</faces>\n");
    fprintf(fp,"</faceSet>\n");
    fprintf(fp,"</geometry>\n");
    fprintf(fp,"</geometries>\n");
    fprintf(fp,"</jvx-model>\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*              FUNCTION TO CALULATE CURVATURE AND OTHER RELATED QUANTITIES                                           */
/*--------------------------------------------------------------------------------------------------------------------*/

void CURVCALC(int num,vertex *vptr,triangle *trptr)
{
    int i,j,im1,ip1,tr1,tr2,k;
    double fwei,nor[4],ssm,pmat[4][4],v2[4],elen,euv[4],fn1[4],fn2[4];
    double an,v1[4],si,xc,dih,en[4],wei,emcu,cpm[3],rtan[3],cmat[4][4];
    double amat[3],smat[3],zdir[3],hmat[4][4],htmat[4][4],chmat[4][4];
    double amsq,smsq,p,q,r,s,insq,e1,e2,evec[4][4],sum;
     zdir[1]=0.00000; zdir[2]=0.00000; zdir[3]=1.00000;
     fwei=0.000; 
     for(i=1;i<=3;i++){
     for(j=1;j<=3;j++)cmat[i][j]=0.000;hmat[i][j]=0.000;pmat[i][j]=0.000;
                      chmat[i][j]=0.000;evec[i][j]=0.0;}

     for(i=1;i<=3;i++){nor[i]=0.0;v2[i]=0.0;euv[i]=0.0;fn1[i]=0.0;fn2[i]=0.0;
                       v1[i]=0.0; en[i]=0.0;cpm[i]=0.0;rtan[i]=0.0;amat[i]=0.0;
                       smat[i]=0.0;}
     weight_calc:    
     for(i=1;i<=((vptr+num)->nonei);i++)                                          /* total area of all nei triangles */
              { j=(vptr+num)->vneitr[i];
                fwei=fwei+(trptr+j)->ar;}
     normal:
     for(i=1;i<=((vptr+num)->nonei);i++){
            j=(vptr+num)->vneitr[i];
            nor[1]=nor[1]+ (((trptr+j)->ar)/fwei)*(trptr+j)->fnor[1];
            nor[2]=nor[2]+ (((trptr+j)->ar)/fwei)*(trptr+j)->fnor[2];
            nor[3]=nor[3]+ (((trptr+j)->ar)/fwei)*(trptr+j)->fnor[3];}
     ssm=sqrt(nor[1]*nor[1]+nor[2]*nor[2]+nor[3]*nor[3]);                       /* normalise the calculated normal   */
     nor[1]=nor[1]/ssm; nor[2]=nor[2]/ssm; nor[3]=nor[3]/ssm;                   /* global frame                      */
     (vptr+num)->vnor[1]=nor[1];                                   
     (vptr+num)->vnor[2]=nor[2];
     (vptr+num)->vnor[3]=nor[3];
     for(i=1;i<=3;i++){
     for(j=1;j<=3;j++){
                  if(i==j) pmat[i][j]=1.0-nor[i]*nor[j];
                  else pmat[i][j]=-1*nor[i]*nor[j];}}
     for(i=1;i<=((vptr+num)->nonei);i++){
             im1=i-1; if(i==1) im1=(vptr+num)->nonei;   
             ip1=i+1; if(i==(vptr+num)->nonei) ip1=1;
                                                            
             v2[1]= ((vptr+((vptr+num)->vneipt[i]))->vcoord[1])-((vptr+num)->vcoord[1]);
             v2[2]= ((vptr+((vptr+num)->vneipt[i]))->vcoord[2])-((vptr+num)->vcoord[2]);
             v2[3]= ((vptr+((vptr+num)->vneipt[i]))->vcoord[3])-((vptr+num)->vcoord[3]);

             elen=sqrt(v2[1]*v2[1]+v2[2]*v2[2]+v2[3]*v2[3]);                   /* length of edge vector from vernu->i */

             v2[1]=v2[1]/elen; v2[2]=v2[2]/elen; v2[3]=v2[3]/elen;             /*  unit vector along the edge        */
             euv[1]=v2[1];euv[2]=v2[2];euv[3]=v2[3];

             tr1=(vptr+num)->vneitr[im1]; tr2=(vptr+num)->vneitr[i];           /*  2 triangles sharing the edge      */

             fn1[1]=(trptr+tr1)->fnor[1]; fn1[2]=(trptr+tr1)->fnor[2]; 
             fn1[3]=(trptr+tr1)->fnor[3];
             fn2[1]=(trptr+tr2)->fnor[1]; fn2[2]=(trptr+tr2)->fnor[2]; 
             fn2[3]=(trptr+tr2)->fnor[3];

             an=0.00000;
             for(j=1;j<=3;j++) an=an+fn1[j]*fn2[j];
             if(an>1.000000000) an=1.000000000000;

             v1[1]=fn1[2]*fn2[3]-fn1[3]*fn2[2];                                /* calculate N1 X N2                 */                
             v1[2]=fn1[3]*fn2[1]-fn1[1]*fn2[3];
             v1[3]=fn1[1]*fn2[2]-fn1[2]*fn2[1];        
     
             si=1.00000;
             xc=0.00000;
             for(j=1;j<=3;j++) xc=xc+euv[j]*v1[j];                             /* E. ( N1 X N2)                     */  
             if(xc<0.00000000000) si=-1.00000;                                 /* sign of E.( N1 X N2 )             */
             dih=pi+si*acos(an);                                               /* signed dihedral angle             */
             en[1]=(trptr+tr1)->fnor[1]+(trptr+tr2)->fnor[1];                  /* normal along above edge           */
             en[2]=(trptr+tr1)->fnor[2]+(trptr+tr2)->fnor[2];
             en[3]=(trptr+tr1)->fnor[3]+(trptr+tr2)->fnor[3];
             ssm=sqrt(en[1]*en[1]+en[2]*en[2]+en[3]*en[3]);
             en[1]=en[1]/ssm;en[2]=en[2]/ssm;en[3]=en[3]/ssm;
             wei=0.0;
             for(j=1;j<=3;j++) wei=wei+nor[j]*en[j];
             emcu=2*elen*cos(dih*0.5);                                         /* mean curvature along the edge     */
 
             cpm[1]=(euv[2]*en[3]-euv[3]*en[2]);                               /* direction orthogonal to edge and  */
             cpm[2]=(euv[3]*en[1]-euv[1]*en[3]);                               /* edge normal                       */
             cpm[3]=(euv[1]*en[2]-euv[2]*en[1]);
   
             for(k=1;k<=3;k++){
             rtan[k]=0.000;
             for(j=1;j<=3;j++) rtan[k]=rtan[k]+pmat[k][j]*cpm[j];}

             ssm=sqrt(rtan[1]*rtan[1]+rtan[2]*rtan[2]+rtan[3]*rtan[3]);
             rtan[1]=rtan[1]/ssm;rtan[2]=rtan[2]/ssm;rtan[3]=rtan[3]/ssm;
             for(k=1;k<=3;k++){
             for(j=1;j<=3;j++)cmat[k][j]=cmat[k][j]+0.5*wei*emcu*rtan[k]*rtan[j];}
        }
        amat[1]=zdir[1]+nor[1];amat[2]=zdir[2]+nor[2];amat[3]=zdir[3]+nor[3];
        smat[1]=zdir[1]-nor[1];smat[2]=zdir[2]-nor[2];smat[3]=zdir[3]-nor[3];
        amsq=amat[1]*amat[1]+amat[2]*amat[2]+amat[3]*amat[3];
        smsq=smat[1]*smat[1]+smat[2]*smat[2]+smat[3]*smat[3];
                     
        house_holder:if(amsq>smsq){                                          /* house holder trans to reduce to 2X2 */
                 amat[1]=amat[1]/sqrt(amsq);amat[2]=amat[2]/sqrt(amsq);amat[3]=amat[3]/sqrt(amsq);
                 for(i=1;i<=3;i++){
                 for(j=1;j<=3;j++){
                 if(i==j) hmat[i][j]=-(1-2*amat[i]*amat[j]);
                 else     hmat[i][j]=2*amat[i]*amat[j];}}}
        else{
                 smat[1]=smat[1]/sqrt(smsq);smat[2]=smat[2]/sqrt(smsq);smat[3]=smat[3]/sqrt(smsq);
                 for(i=1;i<=3;i++){
                 for(j=1;j<=3;j++){
                 if(i==j) hmat[i][j]=(1-2*smat[i]*smat[j]);
                 else     hmat[i][j]=-2*smat[i]*smat[j];}}}
        for(i=1;i<=3;i++){ for(j=1;j<=3;j++)pmat[i][j]=0.0;chmat[i][j]=0.0;}
             
        for(i=1;i<=3;i++){
        for(j=1;j<=3;j++){
        sum=0.000;
        for(k=1;k<=3;k++) {sum=sum+cmat[i][k]*hmat[k][j];                  /* diagonalise the constructed matrix*/
                           chmat[i][j]=sum;}}}

        for(i=1;i<=3;i++){
        for(j=1;j<=3;j++)htmat[i][j]=hmat[j][i];}

        for(i=1;i<=3;i++){
        for(j=1;j<=3;j++){
        sum=0.0;
        for(k=1;k<=3;k++){ sum=sum+htmat[i][k]*chmat[k][j];
                           pmat[i][j]=sum;}}}
        p=pmat[1][1]; q=pmat[1][2]; r=pmat[2][1]; s=pmat[2][2];          /* components of 2 X 2 minor          */
        non_diagonal: if((q !=0.000000) && (r!=0.000000000)){            /* non-diagonal matrices eigen values */
                     insq=(p+s)*(p+s)-4*(p*s-q*r);
                     if ((insq > 0.0)){                                  /* complex values are avoided         */                
                            if((p+s)<0){                                 /* pick up the largest eigen value    */                 
                            e1=((p+s)-sqrt(insq))*0.5;                          
                            e2=((p+s)+sqrt(insq))*0.5;}
                            else{
                            e1=((p+s)+sqrt(insq))*0.5;                   /* eigen values, e1 is the largest    */                    
                            e2=((p+s)-sqrt(insq))*0.5;}}                         
                     else{
                            e1=(p+s)*0.5;                                /* degenerate eigen values            */                 
                            e2=e1;}}
            else{                                                        /* pick up largest eigen value for    */
                 if(fabs(p)>fabs(s)){                                    /* diagonal matrix                    */
                          e1=p ; e2=s;
                          for(i=1;i<=3;i++){for(j=1;j<=3;j++){
                                 if(i==j) evec[i][j]=1;                  /* eigen vector is same as unit matrix*/                   
                                 else evec[i][j]=0;}}}
                else{ e1=s ; e2=p;}}
		 cout << e1 << endl;
		 cout << e2 << endl;
		 cout << (vptr+num)->totarea << endl;
         if(fabs(e1)<pow(10,-10)) e1=0;
         if(fabs(e2)<pow(10,-10)) e2=0;
         (vptr+num)->cur1=e1/((vptr+num)->totarea);                      /* principal curvature 1              */                 
         (vptr+num)->cur2=e2/((vptr+num)->totarea);                      /* principal curvature 2              */
         (vptr+num)->mcur=((vptr+num)->cur1+(vptr+num)->cur2)*0.5;       /* mean curvature                     */
}
/*-------------------------------------------------------------------------------------------------------------------*/
/*                            FUNCTION TO MOVE THE VERTEX                                                            */
/*-------------------------------------------------------------------------------------------------------------------*/
void MOVEVERTEX(int vert,vertex *vptr,triangle *trptr,tether *liptr,tetherm *limptr)
{
   int i,j,k,trno,tn,check;
   double drl,bl,dE;
   double dr[4],co3[4],inen,fien,prob;
   double otar,ocd[4],otr_area[10],otr_vol[10],otr_fnor[10][4], ocu[4]; /* Used for storing the old values  */
   double ovnor[4],ovn_curv[10][4],ovn_vnor[10][4],ovn_totar[10];

   for(i=1;i<=3;i++)  dr[i]=(1.0-2*ran2(&seed))*0.05;  

   ocd[1]=(vptr+vert)->vcoord[1];                                       /* Store old coordinates */
   ocd[2]=(vptr+vert)->vcoord[2];
   ocd[3]=(vptr+vert)->vcoord[3];

   drl=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];
   if(drl>0.0025) return;                                              /* Max allowed displacement is sqrt(0.0025)*/

   (vptr+vert)->vcoord[1]=(vptr+vert)->vcoord[1]+dr[1];                /* New displaced position */
   (vptr+vert)->vcoord[2]=(vptr+vert)->vcoord[2]+dr[2];
   (vptr+vert)->vcoord[3]=(vptr+vert)->vcoord[3]+dr[3];

   nei_len_check:
   for(i=1;i<=((vptr+vert)->nonei);i++){                               /* check minm & maxm distance to neighbhors*/
        co3[1]=(vptr+vert)->vcoord[1]-(vptr+((vptr+vert)->vneipt[i]))->vcoord[1];
        co3[2]=(vptr+vert)->vcoord[2]-(vptr+((vptr+vert)->vneipt[i]))->vcoord[2];
        co3[3]=(vptr+vert)->vcoord[3]-(vptr+((vptr+vert)->vneipt[i]))->vcoord[3];
        bl=sqrt(co3[1]*co3[1]+co3[2]*co3[2]+co3[3]*co3[3]);
        if(bl>sqrt(3.000)||(bl<1.0000000)){ 
                  (vptr+vert)->vcoord[1]=ocd[1];
                  (vptr+vert)->vcoord[2]=ocd[2];
                  (vptr+vert)->vcoord[3]=ocd[3];
                   return;}}

   self_avoid_check:
   for(i=1;i<=nver;i++){                                              /* check minm distance to neighbhors*/
   if(i!=vert){ co3[1]=(vptr+vert)->vcoord[1]-((vptr+i)->vcoord[1]);
                 co3[2]=(vptr+vert)->vcoord[2]-((vptr+i)->vcoord[2]);
                 co3[3]=(vptr+vert)->vcoord[3]-((vptr+i)->vcoord[3]);
                 bl=sqrt(co3[1]*co3[1]+co3[2]*co3[2]+co3[3]*co3[3]);
                 if(bl<1.0000000){
                      (vptr+vert)->vcoord[1]=ocd[1];
                      (vptr+vert)->vcoord[2]=ocd[2];
                      (vptr+vert)->vcoord[3]=ocd[3];
                      return;}}}

    otar=(vptr+vert)->totarea;                                        /* save old values of vertex & neighbhors*/
    nei_ar: for(i=1;i<=((vptr+vert)->nonei);i++){ 
                   j=(vptr+vert)->vneipt[i];                          /* neighbouring vertex */
                   k=(vptr+vert)->vneitr[i];                          /* neighbouring triangle */                        
                   otr_area[i]=(trptr+k)->ar;                         /* otr - old triangle */ 
                   otr_vol[i]=(trptr+k)->vol;                         /* ovn- old vertex*/
                   ovn_totar[i]=(vptr+j)->totarea; 
                   otr_fnor[i][1]=(trptr+k)->fnor[1];
                   otr_fnor[i][2]=(trptr+k)->fnor[2];
                   otr_fnor[i][3]=(trptr+k)->fnor[3];             
                   }
      inen=0.0;
      mv_neigh: for (i=1;i<=((vptr+vert)->nonei);i++){                /* initial energy calculation*/
                   j=(vptr+vert)->vneipt[i];                          /* coordinate changes does not matter here  */
                   k=(vptr+vert)->vneitr[i];                          /* so far only coord of vert changed to new */
                   inen=inen+(((vptr+j)->mcur)*((vptr+j)->mcur))*
                   kappa*((vptr+j)->totarea)-pr*(trptr+k)->vol;}      
      inen=inen+(((vptr+vert)->mcur)*((vptr+vert)->mcur))*kappa
                 *((vptr+vert)->totarea);

      upd_ntr_area:
      for(i=1;i<=((vptr+vert)->nonei);i++) 
      areacal((vptr+vert)->vneitr[i],trptr,vptr);                     /* new area of nei vertex and area,volume &*/
                                                                      /* normal of nei triangle calculated here  */        
      face_ang_check:check=0; 
      for(trno=1;trno<=((vptr+vert)->nonei);trno++){                  /* select one nei- triangle                */          
            tn=(vptr+vert)->vneitr[trno];                         
            check=faceangchk(tn,trptr,liptr,limptr);
            if (check!=0){                                            /* maxm allowed angle is 0.5 */        
                         (vptr+vert)->vcoord[1]=ocd[1];     
                         (vptr+vert)->vcoord[2]=ocd[2]; 
                         (vptr+vert)->vcoord[3]=ocd[3];
                         (vptr+vert)->totarea=otar;
                         for(i=1;i<=((vptr+vert)->nonei);i++){
                                j=(vptr+vert)->vneipt[i];  
                 		k=(vptr+vert)->vneitr[i]; 
   				(trptr+k)->ar=otr_area[i];
   		                (vptr+j)->totarea=ovn_totar[i];
                                (trptr+k)->fnor[1]=otr_fnor[i][1];
                                (trptr+k)->fnor[2]=otr_fnor[i][2];
                                (trptr+k)->fnor[3]=otr_fnor[i][3];}    /* restore volume here                   */
                                return; 
        }}
 
     ocu[1]=(vptr+vert)->cur1;  ovnor[1]=(vptr+vert)->vnor[1];         /* store current properties of moved     */
     ocu[2]=(vptr+vert)->cur1;  ovnor[2]=(vptr+vert)->vnor[2];         /* vertrx and its neighbours             */
     ocu[3]=(vptr+vert)->mcur;  ovnor[3]=(vptr+vert)->vnor[3];              
     for(i=1;i<=((vptr+vert)->nonei);i++){
           j=(vptr+vert)->vneipt[i];
           ovn_curv[i][1]=((vptr+j)->cur1);
           ovn_curv[i][2]=((vptr+j)->cur2);
           ovn_curv[i][3]=((vptr+j)->mcur);
           ovn_vnor[i][1]=(vptr+j)->vnor[1];
           ovn_vnor[i][2]=(vptr+j)->vnor[2];
           ovn_vnor[i][3]=(vptr+j)->vnor[3];}

     CURVCALC(vert,vptr,trptr);                                        /*curvature calculation of vertex and nei*/
     for(i=1;i<=((vptr+vert)->nonei);i++)
     CURVCALC( (vptr+vert)->vneipt[i],vptr,trptr);

     fien=0.0000;                                                      /* final energy calculation                */
     for (i=1;i<=((vptr+vert)->nonei);i++){
                   j=(vptr+vert)->vneipt[i]; 
                   k=(vptr+vert)->vneitr[i];
                   fien=fien+(((vptr+j)->mcur)*((vptr+j)->mcur))
                   *kappa*((vptr+j)->totarea)-pr*(trptr+k)->vol;}
     fien=fien+(((vptr+vert)->mcur)*((vptr+vert)->mcur))*kappa*
               ((vptr+vert)->totarea);

     dE=fien-inen;   
     mvmetro1: if(dE>0.0){                                              /* metropolis                             */
     prob=exp(-beta*dE);
     mvmetro2: if(prob<ran2(&seed)){                                    /* if metropolis fails                    */                              (vptr+vert)->vcoord[1]=ocd[1];                             /* restore old properties of moved vertex */
             (vptr+vert)->vcoord[2]=ocd[2];
             (vptr+vert)->vcoord[3]=ocd[3];
             (vptr+vert)->totarea=otar;
             (vptr+vert)->cur1=ocu[1];
             (vptr+vert)->cur2=ocu[2];
             (vptr+vert)->mcur=ocu[3];
             (vptr+vert)->vnor[1]=ovnor[1];
             (vptr+vert)->vnor[2]=ovnor[2];
             (vptr+vert)->vnor[3]=ovnor[3];
             for (i=1;i<=((vptr+vert)->nonei);i++){                     /* restore old properties of nei          */
                      j=(vptr+vert)->vneipt[i]; 
                      k=(vptr+vert)->vneitr[i];
 
                      (trptr+k)->ar=otr_area[i];
                      (trptr+k)->vol=otr_vol[i];

   		      (vptr+j)->totarea=ovn_totar[i];
                      (vptr+j)->cur1=ovn_curv[i][1];
                      (vptr+j)->cur2=ovn_curv[i][2];
                      (vptr+j)->mcur=ovn_curv[i][3];

                      (trptr+k)->fnor[1]=otr_fnor[i][1];
                      (trptr+k)->fnor[2]=otr_fnor[i][2];
                      (trptr+k)->fnor[3]=otr_fnor[i][3];

                      (vptr+j)->vnor[1]=ovn_vnor[i][1];
                      (vptr+j)->vnor[2]=ovn_vnor[i][2];
                      (vptr+j)->vnor[3]=ovn_vnor[i][3];}}
  }
}
/*-------------------------------------------------------------------------------------------------------------------------*/
/*                                        FUNCTION TO FLIP A tether                                                          */
/*-------------------------------------------------------------------------------------------------------------------------*/

void FLIPPING(int rand,vertex *vptr,triangle *trptr,tether *liptr,tetherm *limptr)
{

int i,j,k,m1,m2,lp1,lp2,tn1,tn2,vt1,vt2,ep1,ep2,fv1,fv2,v,ch1,check,mchk;
int t1[4],t2[4],cver[5],tmp1[4],tmp[4],fvp1,fvm1,fvp2,fvm2,trn;
int tdt_vert[3][4],tdt_li[3][4],tdv_nonnei[5],tdv_neitr[5][11],tdv_neipt[5][11];
double tdt_ar[3],tdt_vol[3],tdt_fnor[3][4],tdv_cur1[5],tdv_cur2[5],tdv_mcur[5];
double tdv_vnor[5][4],tdv_totarea[5];
double inen,fien,dE,prob,rnum,bl,co3[4];


       if(rand>0){                                                           /* t1 & t2 are vertices tethered to the */
           for(i=1;i<=3;i++){ t1[i]=(trptr+((liptr+rand)->tr))->vert[i];     /* triangle on both sides of tether     */
                              t2[i]=(trptr+((limptr+rand)->tr))->vert[i];}
           tn1=(liptr+rand)->tr; tn2=(limptr+rand)->tr;                      /* triangle on both sides of tether     */
           lp1=(liptr+rand)->sep[1]; lp2=(liptr+rand)->sep[2];}              /* end points of tether                 */
       else{
           for(i=1;i<=3;i++){ t1[i]=(trptr+((limptr-rand)->tr))->vert[i];
                              t2[i]=(trptr+((liptr-rand)->tr))->vert[i];}

           t2[3]=(trptr+((liptr-rand)->tr))->vert[3];
           tn1=(limptr-rand)->tr; tn2=(liptr-rand)->tr;
           lp1=(limptr-rand)->sep[1]; lp2=(limptr-rand)->sep[2];}
  
      ep1=0;ep2=0;
      for(i=1;i<=3;i++){
           if (t1[i]==lp2) ep1=i;
           if (t2[i]==lp1) ep2=i;
         
           if((t1[i]!=lp1)&& (t1[i]!=lp2)){ vt1=t1[i]; fv1=i; }             /* lp1 and lp2 are tethered points and     */                 if((t2[i]!=lp2)&&(t2[i]!=lp1)){ vt2=t2[i]; fv2=i;} }                     /* vt1 and vt2 are free vertices         */                                                                                          /* fv1 and fv2 are the position of the   */
                                                                            /* free vertices                         */
      cver[1]=vt1; cver[2]=vt2; cver[3]=lp1; cver[4]=lp2;                   /* All relevant vertice num in one array */
     for(i=1;i<=3;i++){
       tdt_vert[1][i]=(trptr+tn1)->vert[i];                                 /* To store the initial attributes       */
       tdt_vert[2][i]=(trptr+tn2)->vert[i];
       tdt_li[1][i]=(trptr+tn1)->li[i];
       tdt_li[2][i]=(trptr+tn2)->li[i];
       tdt_fnor[1][i]=(trptr+tn1)->fnor[i];
       tdt_fnor[2][i]=(trptr+tn2)->fnor[i];}
       tdt_ar[1]=(trptr+tn1)->ar;
       tdt_ar[2]=(trptr+tn2)->ar;
       tdt_vol[1]=(trptr+tn1)->vol;
       tdt_vol[2]=(trptr+tn2)->vol;

    for(i=1;i<=4;i++){
              tdv_nonnei[i]=(vptr+cver[i])->nonei;
              for(j=1;j<=10;j++){tdv_neitr[i][j]=(vptr+cver[i])->vneitr[j];
                                 tdv_neipt[i][j]=(vptr+cver[i])->vneipt[j];}
               tdv_cur1[i]=(vptr+cver[i])->cur1;
               tdv_cur2[i]=(vptr+cver[i])->cur2;
               tdv_mcur[i]=(vptr+cver[i])->mcur;
               tdv_totarea[i]=(vptr+cver[i])->totarea;
               for(j=1;j<=3;j++) tdv_vnor[i][j]=(vptr+cver[i])->vnor[j];}

       for(i=1;i<=((vptr+vt1)->nonei);i++){
                if (((vptr+vt1)->vneipt[i])==vt2) return;}                    /* Do not proceed if already connected */
       bl_check:
       co3[1]=(vptr+vt2)->vcoord[1]-((vptr+vt1)->vcoord[1]);                  /* check for new bond length          */
       co3[2]=(vptr+vt2)->vcoord[2]-((vptr+vt1)->vcoord[2]);
       co3[3]=(vptr+vt2)->vcoord[3]-((vptr+vt1)->vcoord[3]);
       bl=sqrt(co3[1]*co3[1]+co3[2]*co3[2]+co3[3]*co3[3]);
       if((bl>sqrt(3.00000))||(bl<=1.0000))  return;

       if((((vptr+vt1)->nonei)>=9)||(((vptr+vt2)->nonei)>=9))return;          /* Maxm & minm limit on no.of neighbors*/
       if((((vptr+lp1)->nonei)<=3)&&(((vptr+lp2)->nonei)<=3))return;      

       inen=0.0;
       four_points:for(i=1;i<=4;i++){                                        /* initial energy calculation */
           v=cver[i]; 
           inen=inen+(((vptr+v)->mcur)*((vptr+v)->mcur))
                      *kappa*((vptr+v)->totarea);}
       inen=inen-pr*(((trptr+tn1)->vol)+((trptr+tn2)->vol)) ;
 
       t1[ep1]=vt2; t2[ep2]=vt1;                                             /* free vertices are connected now */
       j=(vptr+cver[3])->nonei;

       for(i=1;i<=j-1;i++){                                                  /* lp2 is removed from the list of lp1 */
       if (((vptr+cver[3])->vneipt[i])==cver[4]){
       for(k=i;k<=j-1;k++) (vptr+cver[3])->vneipt[k]=
                            (vptr+cver[3])->vneipt[k+1];}}
       for(k=j;k<=10;k++)((vptr+cver[3])->vneipt[k])=0;                      /* all sites beyond neighbour size = 0 */
       (vptr+cver[3])->nonei=j-1;
       j=(vptr+cver[4])->nonei;

       for(i=1;i<=j-1;i++){                                                  /* lp1 is removed from the list of lp2 */
       if( ((vptr+cver[4])->vneipt[i])==cver[3]){
       for(k=i;k<=j-1;k++) (vptr+cver[4])->vneipt[k]=
                       (vptr+cver[4])->vneipt[k+1];}}
       for(k=j;k<=10;k++)((vptr+cver[4])->vneipt[k])=0;                      /* all sites beyond neighbour size = 0 */
       (vptr+cver[4])->nonei=j-1;

       tmp1[1]=0;tmp1[2]=0;tmp1[3]=0;
       for(i=1;i<=3;i++){                                                    /* To find the tethers whose triangles   */ 
          m1= (trptr+tn1)->li[i];                                            /* will change                         */
          if(m1>0){if ((liptr+m1)->sep[2]==vt1) tmp1[1]=(trptr+tn1)->li[i];} /* update start and end of chosen tether */
          else {if((limptr+(-m1))->sep[2]==vt1) tmp1[1]=(trptr+tn1)->li[i];}

          m2=(trptr+tn2)->li[i];
          if(m2>0){ if ((liptr+m2)->sep[2]==vt2) tmp1[2]=(trptr+tn2)->li[i];}
          else {if((limptr+(-m2))->sep[2]==vt2) tmp1[2]=(trptr+tn2)->li[i];}}

       if(rand>0) { (liptr+rand)->sep[1]=vt2; (liptr+rand)->sep[2]=vt1;
                    (limptr+rand)->sep[1]=vt1; (limptr+rand)->sep[2]=vt2;}

       else  { (limptr+(-rand))->sep[1]=vt2; (limptr+(-rand))->sep[2]=vt1;
               (liptr+(-rand))->sep[1]=vt1; (liptr+(-rand))->sep[2]=vt2;}

       for (i=1;i<=3;i++){(trptr+tn1)->vert[i]=t1[i];                        /* update the vertex of triangles     */
                          (trptr+tn2)->vert[i]=t2[i];}

       fvp1=fv1+1 ; if(fv1==3)fvp1=1;                                        /* circular boundary conditions       */ 
       fvm1=fv1-1 ; if(fv1==1)fvm1=3;
       fvp2=fv2+1 ; if(fv2==3)fvp2=1;
       fvm2=fv2-1 ; if(fv2==1)fvm2=3;

       tmp[1]=(trptr+tn1)->li[fvm1];                                         /* tmp is a temp array to stores tethers*/
       tmp[2]=(trptr+tn2)->li[fvm2];

       (trptr+tn1)->li[fvm1]=(trptr+tn1)->li[fvp1];                          /* updating the triangle -- tether      */ 
       (trptr+tn2)->li[fvm2]=(trptr+tn2)->li[fvp2];
       (trptr+tn1)->li[fvp1]=tmp[2];
       (trptr+tn2)->li[fvp2]=tmp[1];

       
       if(tmp1[1]>0) (liptr+tmp1[1])->tr=tn2 ;                              /* updating tether -- triangle            */
       else (limptr+(-tmp1[1]))->tr=tn2 ;
       if(tmp1[2]>0) (liptr+tmp1[2])->tr=tn1 ;
       else (limptr+(-tmp1[2]))->tr=tn1 ;

       (vptr+cver[1])->nonei=(vptr+cver[1])->nonei+1;                       /* increment of neigh of ver1&ver2 by 1 */
       (vptr+cver[2])->nonei=(vptr+cver[2])->nonei+1;

       ch1=0 ; i=1;
       f1_ver:while((i<=tdv_nonnei[1])&& (ch1==0)){                         /* include change in nei of vertex      */
              if(tdv_neipt[1][i] == cver[3]){
               for(j=i+1;j<=tdv_nonnei[1];j++){
              (vptr+cver[1])->vneipt[j+1]=tdv_neipt[1][j];}                    
              (vptr+cver[1])->vneipt[i+1]=cver[2] ; ch1=1;}
               i=i+1;}

         ch1=0 ; i=1;
         f2_ver:while((i<=tdv_nonnei[2])&& (ch1==0)){        
              if(tdv_neipt[2][i]==cver[4]){
              for(j=i+1;j<=tdv_nonnei[2];j++)
               (vptr+cver[2])->vneipt[j+1]=tdv_neipt[2][j];
              (vptr+cver[2])->vneipt[i+1]=cver[1] ;ch1=1;}                        
               i=i+1;}  

         ch1=0 ; i=1;                                                                  
         f1_tr:while((i<=tdv_nonnei[1]) && (ch1==0)){                      /* include change in nei of triangle */    
             if (tdv_neitr[1][i]==tn1){                                              
             for(j=i+1;j<=tdv_nonnei[1];j++) 
             (vptr+cver[1])->vneitr[j+1]=tdv_neitr[1][j];
             (vptr+cver[1])->vneitr[i+1]=tn2 ;ch1=1;}
             i=i+1;}

         ch1=0 ; i=1;
         f2_tr: while((i<=tdv_nonnei[2]) && (ch1 ==0)){ 
              if(tdv_neitr[2][i]==tn2){
              for(j=i+1;j<=tdv_nonnei[2];j++) 
              (vptr+cver[2])->vneitr[j+1]=tdv_neitr[2][j];
              (vptr+cver[2])->vneitr[i+1]=tn1 ;ch1=1;}
              i=i+1;}

         ch1=0 ; i=1;
         f3_tr:while(i<=tdv_nonnei[3] && (ch1==0)){          
             if((vptr+cver[3])->vneitr[i]== tn2){                       
             for(j=i;j<=((vptr+cver[3])->nonei);j++)
              (vptr+cver[3])->vneitr[j]=(vptr+cver[3])->vneitr[j+1]; 
             if(i==tdv_nonnei[3]){                                
              (vptr+cver[3])->vneitr[tdv_nonnei[3]]=(vptr+cver[3])->vneitr[1];  
              for(k=1;k<=((vptr+cver[3])->nonei);k++) 
              (vptr+cver[3])->vneitr[k]= (vptr+cver[3])->vneitr[k+1];}
              ch1=1;}
              i=i+1;}
         for(i=tdv_nonnei[3];i<=10;i++)  (vptr+cver[3])->vneitr[i]=0 ;

         ch1=0 ; i=1;
            f4_tr:while((i<=tdv_nonnei[4])&& (ch1==0)){          
             if ((vptr+cver[4])->vneitr[i] == tn1){
             for(j=i;j<=((vptr+cver[4])->nonei);j++) 
             (vptr+cver[4])->vneitr[j]=(vptr+cver[4])->vneitr[j+1];
             if(i == tdv_nonnei[4]){
             (vptr+cver[4])->vneitr[tdv_nonnei[4]]=(vptr+cver[4])->vneitr[1];
             for(k=1;k<=((vptr+cver[4])->nonei);k++)
             (vptr+cver[4])->vneitr[k]= (vptr+cver[4])->vneitr[k+1];
             } 
             ch1=1;
             }
             i=i+1;}
         for(k=tdv_nonnei[4];k<=10;k++) (vptr+cver[4])->vneitr[k]=0;
         onlyarea(tn1,trptr,vptr);                                          /* calculate area of new triangles  */  
         onlyarea(tn2,trptr,vptr);

         for(i=1;i<=4;i++){
             (vptr+cver[i])->totarea=0;
             for( j=1;j<=((vptr+cver[i])->nonei);j++){
             trn=((vptr+cver[i])->vneitr[j]);
             (vptr+cver[i])->totarea=((vptr+cver[i])->totarea)+
                                      ((trptr+trn)->ar)/3.0;}}

         mchk=0;check=0; fien=0.0;
         check=faceangchk(tn1,trptr,liptr,limptr);                       /* check the angle between tri 1 and   */
                                                                         /* neigh faces                         */
         if(check==0){
            check=faceangchk(tn2,trptr,liptr,limptr);
            if(check==0){        
                  for(i=1;i<=4;i++) CURVCALC(cver[i],vptr,trptr);        /* update curvature at chosen vertices */
          
                  for(i=1;i<=4;i++){ v=cver[i];                           
                  fien=fien+(((vptr+v)->mcur)*((vptr+v)->mcur))*         /* final energy calculations           */
                              kappa*((vptr+v)->totarea);}
                  fien=fien-pr*(((trptr+tn1)->vol)+((trptr+tn2)->vol)) ;
                  dE=fien-inen;
                  if(dE>0.00){ prob=exp(-beta*dE);                       /* probability is the Boltzman factor  */                                      rnum=ran2(&seed);                                                                    
                               fmetro2:if(prob<rnum) mchk=1;                         
                                  else mchk=0;}
                            }
              else mchk=1;}
    else mchk=1;
 
    flip_failed: if(mchk==1){
         if(rand>0){                                                     /* Restore the old properties */
         (liptr+rand)->sep[1]=cver[3]; (liptr+rand)->sep[2]=cver[4];
         (limptr+rand)->sep[1]=cver[4] ; (limptr+rand)->sep[2]=cver[3];}
         else{ 
         (limptr-rand)->sep[1]=cver[3]  ; (limptr-rand)->sep[2]=cver[4];         
         (liptr-rand)->sep[1]=cver[4] ; (liptr-rand)->sep[2]=cver[3];}

         for(i=1;i<=3;i++){
             (trptr+tn1)->vert[i]=tdt_vert[1][i];
             (trptr+tn2)->vert[i]=tdt_vert[2][i];             
             (trptr+tn1)->li[i]=tdt_li[1][i]    ;
             (trptr+tn2)->li[i]=tdt_li[2][i];
             (trptr+tn1)->fnor[i]=tdt_fnor[1][i];
             (trptr+tn2)->fnor[i]=tdt_fnor[2][i];}                                
             (trptr+tn1)->ar=tdt_ar[1]     ; (trptr+tn2)->ar=tdt_ar[2];                
             (trptr+tn1)->vol=tdt_vol[1]   ; (trptr+tn2)->vol=tdt_vol[2];

        for(i=1;i<=4;i++){              
            (vptr+cver[i])->nonei=tdv_nonnei[i];                                                             
            (vptr+cver[i])->cur1=tdv_cur1[i];                                              
            (vptr+cver[i])->cur2=tdv_cur2[i];                                               
            (vptr+cver[i])->mcur=tdv_mcur[i];                                       
            (vptr+cver[i])->totarea=tdv_totarea[i]; 
             for(k=1;k<=3;k++) (vptr+cver[i])->vnor[k]=tdv_vnor[i][k];
             for(k=1;k<=9;k++){ 
               (vptr+cver[i])->vneitr[k]=tdv_neitr[i][k];                               
               (vptr+cver[i])->vneipt[k]=tdv_neipt[i][k];}}

        for(i=1;i<=3;i++){
            m1=(trptr+tn1)->li[i] ; m2=(trptr+tn2)->li[i];
            if(m1>0) (liptr+m1)->tr=tn1;
            else  (limptr+(-m1))->tr=tn1;
                       
            if(m2>0)(liptr+m2)->tr=tn2;
            else (limptr+(-m2))->tr=tn2;}
        }                                            
}
/*------------------------------------------------------------------------------------------------------------------*/
/*                                    FUNCTION TO FIND ONLY AREA                                                    */
/*------------------------------------------------------------------------------------------------------------------*/
void onlyarea(int tr,triangle *trptr,vertex *vptr)
{
     int i,j,k,n;
     double r1[4],r2[4],r3[4],r21[4],r31[4],ax,ay,az,area;
     ax=0.0; ay=0.0; az=0.0;
 
     i=(trptr+tr)->vert[1];j=(trptr+tr)->vert[2];k=(trptr+tr)->vert[3];               /* vertices of the triangle   */ 
     for(n=1;n<=3;n++){ r1[n]=(vptr+i)->vcoord[n];                                    /* coordinates of vertices    */
                        r2[n]=(vptr+j)->vcoord[n]; 
                        r3[n]=(vptr+k)->vcoord[n]; }   

     for(n=1;n<=3;n++){ r21[n]=r2[n]-r1[n]; r31[n]=r3[n]-r1[n];} 
     ax=(r21[2]*r31[3]-r21[3]*r31[2])*0.5;
     ay=(r21[3]*r31[1]-r21[1]*r31[3])*0.5;
     az=(r21[1]*r31[2]-r21[2]*r31[1])*0.5;
    
     area=sqrt(ax*ax+ay*ay+az*az);
     (trptr+tr)->ar=area; 
     ax=ax/area ; ay=ay/area ;  az=az/area;
     (trptr+tr)->fnor[1]=ax; (trptr+tr)->fnor[2]=ay; 
     (trptr+tr)->fnor[3]=az;
     (trptr+tr)->vol=(r1[1]*(r2[2]*r3[3]-r2[3]*r3[2])+
                      r1[2]*(r2[3]*r3[1]-r2[1]*r3[3])+
                  r1[3]*(r2[1]*r3[2]-r2[2]*r3[1]))/6.0;
  return;
}
/*------------------------------------------------------------------------------------------------------------------*/
/*                                  FUNCTION FOR FACE ANGLE CHECK                                                   */
/*------------------------------------------------------------------------------------------------------------------*/
int faceangchk(int trno,triangle *trptr,tether *liptr, tetherm *limptr)
{
   int i,j,k,ll[4],check;
   double npro;
   check=0;

   for(i=1;i<=3;i++)ll[i]=(trptr+trno)->li[i];                            /* store it tethers                         */
   for(k=1;k<=3;k++){
          if(ll[k]<0) j=(liptr+(-ll[k]))->tr;                             /* select one of the tethers and find other */
          else j=(limptr+ll[k])->tr;                                      /* triangle connected to it               */
          npro=0.0;
          for(i=1;i<=3;i++)                                               /* find angle between both triangle       */
          npro=npro+((trptr+trno)->fnor[i])*((trptr+j)->fnor[i]); 
          if (npro<fangle) check=1;}
   return(check);
}
/*------------------------------------------------------------------------------------------------------------------*/
/*                                  FUNCTION TO CALCULATE SYSTEM ENERGY                                             */
/*------------------------------------------------------------------------------------------------------------------*/
void system_energy(int n,triangle *trptr,vertex *vptr)
{
  int i,j;
  float elen,volen,tmcur;
  elen=0; volen=0; tmcur=0;

  for(i=1;i<=nver;i++)
  elen=elen+((vptr+i)->mcur)*((vptr+i)->mcur)*kappa*((vptr+i)->totarea);
  
  for(i=1;i<=nver;i++) tmcur=tmcur+(vptr+i)->mcur;

  for(i=1;i<=ntr;i++) volen=volen+(trptr+i)->vol; 
   
  fprintf(fp1,"%d\t%f\t%f\t%f\n",n,elen,volen,tmcur);  
}
/*-----------------------------------------------------------------------------------------------------------------*/


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
             

