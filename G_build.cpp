/*
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 1.
this C++ code is to generate a Resch topological structure, that one 'hex' is degree one,
and extend to seven is degree two.
the input can only be modified by the constant 'Degree'
and the output is 
V,S,E
V-th lines following, each item is:
index of the vertex;
the type of a vertex, that -1 is in the surrounding and others may be useless in a more general usage;
the x position in a unfolded 2D state;
the y position in a unfolded 2D state;
E-th lines following each item is:
index of the edge;
the type of a vertex, that 1 and 2 means different folding protrusion directions and -1 in the surrounding;
the u vertex in the node;
the v vertex in the node, insure that a non-directional edge will only appear once;
the first surface adjacenced by the edge;
the second surface adjacenced by the edge, this can be used to build a Dual Graph of this topology.
NOT a python code is just for that it is completed earlier than others and I have no realization that it can be used in
a larger project as the first step. MEANWHILE C++ is easier than others when pure algorithm is concentrated. 
as a early design, the output file is outGraph.out in '../out/'.
*/
#include<cstdio>
#include<cmath>
#include<iostream>
using namespace std;
#define DEGREE 1 //degree of the topology, of which 0 is the original point
#define VERTEX_N 2000
#define EDGE_N 12000
#define HALFX 0.8660254037844386     //sqrt3 /2
#define HALFY 0.5
#define PI 3.1415926535
struct Edge{
    int u,v; //from u to v
    int id;
    int type;//type1=valley type2=mount
    int split;//if type 2 split e=(u,v), split=v
    int sur=-1;
    Edge *next;
}E[EDGE_N]; int emp;//the last menory of E used
struct Vertex{
    int id;
    int type;
    double x,y;
    int in_deg;
    Edge *head;
}V[VERTEX_N];int vmp;//the last menory of V used
int numS=0;
namespace Gbuild{

int edge_is_existed[6];
bool equal_float(double a,double b){
    #define eps 1e-4
    return abs(a-b)<eps;
}
const double ID_TO_YX[6][2]={1,0,HALFY,-HALFX,-HALFY,-HALFX,-1,0,-HALFY,HALFX,HALFY,HALFX};
int get_rel_id(Vertex *u,Vertex *v){
    double x=v->x-u->x;
    double y=v->y-u->y;
    for(int i=0;i<6;i++){
        if(equal_float(y,ID_TO_YX[i][0])&&equal_float(x,ID_TO_YX[i][1]))return i;
    }
    return -1;
}


Vertex *queue[VERTEX_N];int head,tail;
void pushq(Vertex *v){
    queue[head++]=v;
}
Vertex *popq(){
    return queue[tail++];
}
void newEdge(Vertex *u,Vertex *v){
    emp++;
    E[emp].u=u->id;
    E[emp].v=v->id;
    E[emp].next=u->head;
    E[emp].split=-1;
    E[emp].type=-1;
    u->head=&E[emp];
    E[emp].id=emp;
}
Vertex* newVertex(Vertex *A,int relid){
    vmp++;
    V[vmp].in_deg=A->in_deg+1;
    V[vmp].id=vmp;
    V[vmp].type=1;
    V[vmp].y=A->y+ID_TO_YX[relid][0];
    V[vmp].x=A->x+ID_TO_YX[relid][1];
    pushq(&V[vmp]);
    return &V[vmp];
}
void buildHex(Vertex *A){
    for(int i=0;i<6;i++)edge_is_existed[i]=-1;
    Vertex *Vx;
    for(Edge *cp=A->head;cp;cp=cp->next){
        edge_is_existed[get_rel_id(A,&V[cp->v])]=cp->v;
    }
    for(int i=0;i<6;i++){
        if(edge_is_existed[i]==-1){
            Vx=newVertex(A,i);
            edge_is_existed[i]=Vx->id;
            newEdge(A,Vx);
            newEdge(Vx,A);
        }
    }
    bool isadj;
    for(int i=0;i<6;i++){
        isadj=false;
        for(Edge *cp=V[edge_is_existed[i]].head;cp;cp=cp->next){
            if(cp->v==edge_is_existed[(i+1)%6])isadj=true;
        }
        if(!isadj){
            newEdge(&V[edge_is_existed[i]],&V[edge_is_existed[(i+1)%6]]);
            newEdge(&V[edge_is_existed[(i+1)%6]],&V[edge_is_existed[i]]);
        }
    }
}

void buildHex0(){
    V[0].type=1;//all is 0 except type
    pushq(&V[0]);
    Vertex *A;
    for(int i=0;i<DEGREE;i++){
        while(queue[tail]->in_deg==i){
            A=popq();
            buildHex(A);
        }
    }
}
//attention: split item of split vertex is meaningless
void newEdge(Vertex *u,Vertex *v,int to_split,int type){
    newEdge(u,v);
    E[emp].split=to_split;
    E[emp].type=type;
}

Vertex *split(Vertex *u,Vertex *v,int t){
    vmp++;
    V[vmp].id=vmp;
    V[vmp].type=2;
    V[vmp].y=(v->y-u->y)*(t^3)/3+u->y;
    V[vmp].x=(v->x-u->x)*(t^3)/3+u->x;
    if(u->in_deg==DEGREE&&v->in_deg==DEGREE)V[vmp].type=-1;
    newEdge(u,&V[vmp],v->id,t);
    newEdge(&V[vmp],u,-1,t);
    newEdge(v,&V[vmp],u->id,t^3);
    newEdge(&V[vmp],v,-1,t^3);//if t=2 t'=1 or otherwise
    return &V[vmp];
}


Vertex *vertexadj[6];
void getSplit(Vertex *A){
    for(int i=0;i<6;i++)vertexadj[i]=0;
    int relid;
    for(Edge *cp=A->head;cp;cp=cp->next){
        relid=get_rel_id(A,&V[cp->v]);
        if(relid==-1){//relid =-1 means that this is a directorical vertex split
            relid=get_rel_id(A,&V[cp->split]); 
            vertexadj[relid]=&V[cp->v];
        }else{
            if(!vertexadj[relid])vertexadj[relid]=&V[cp->v];
        }
        
    }
    int etype0;//at 1/3(type2) or 2/3(type1)
    for(int i=0;i<6;i++){
        if(vertexadj[i]){
            if(vertexadj[i]->type==1){
                etype0=i%2+1;
                vertexadj[i]=split(A,vertexadj[i],etype0);
            }
        }
    }
    for(int i=0;i<6;i++){
        if(vertexadj[i]&&vertexadj[(i+1)%6]){
            newEdge(vertexadj[i],vertexadj[(i+1)%6],-1,2);
            newEdge(vertexadj[(i+1)%6],vertexadj[i],-1,2);
        }
    }
}
void getSplit0(){
    int vmp0=vmp;
    for(int i=0;i<=vmp0;i++){
        getSplit(&V[i]);
    }
}
}
namespace Dual{
int visited[EDGE_N];

void wrap(int n,Edge *p){
    p->sur=n;
    visited[p->id]=1;
    double theta=atan2(V[p->u].y-V[p->v].y,V[p->u].x-V[p->v].x),eta;
    double min=1e8;
    Edge *left=0;
    for(Edge *q=V[p->v].head;q;q=q->next){
        if(visited[q->id]==-1)continue;
        eta=atan2(V[q->v].y-V[q->u].y,V[q->v].x-V[q->u].x);
        eta=theta-eta;
        if(eta<0)eta=2*PI+eta;
        if(q->v!=p->u&&eta<min){
            min=eta;
            left=q;
        }
    }
    if(visited[left->id]==0)wrap(n,left);

}
void getsurface(){
    for(int i=0;i<=emp;i++){
        if(E[i].type==-1)visited[i]=-1;
    }
    for(int i=1;i<=emp;i++){
        if(visited[i]==0){    
            wrap(numS,&E[i]);
            numS++;
        }
    }
    //here to tick the most outscale , for it's unexisted, sur=-1
    for(int i=0;i<numS;i++)visited[i]=0;
    for(int i=1;i<=emp;i++){
        if(E[i].type==-1)continue;
        visited[E[i].sur]++;
    }
    for(int i=1;i<=emp;i+=2){
        if(visited[E[i].sur]>3){
            E[i].sur=E[i+1].sur;
            E[i+1].sur=-1;
        }
        if(visited[E[i+1].sur]>3){
            E[i+1].sur=-1;
        }
    }
    //we now found that use this tratics, the outter scale will only be the last one, so exclude it is samply numS--
    numS-=1;
    
}

} 

int main(){
    printf("if want to change the degree, do it in #define, and fix the memory size of V and E.\n");
    printf("the output is the Graph of the topology of the structure, where E is non-dir.\n");
    printf("first line of output is num of V, S and E, and the next n_V lines are Vi's id/type/x/y, and the next n_E lines are Ei's id/type/u/v/nearS1/nearS2.\n");

    Gbuild::buildHex0();
    printf("%d %d\n",vmp,emp);
    int emp0=emp;//where emp0 is the num of E unsplit.
    Gbuild::getSplit0();
    printf("%d %d",vmp,emp);
    Dual::getsurface();
    FILE *out=fopen("out\\outGraph.out","w");
    fprintf(out,"%d %d %d\n",vmp+1,numS,(emp-emp0)/2);
    for(int i=0;i<=vmp;i++){
        if(V[i].in_deg==DEGREE)V[i].type=-1;
        fprintf(out,"%d %d %.12lf %.12lf\n",V[i].id,V[i].type,V[i].x,V[i].y);
    }
    //emp is always pairly created, so we cut half.
    for(int i=1;i<=emp;i+=2){
        if(E[i].type!=-1){
            fprintf(out,"%d %d %d %d %d %d\n",(E[i].id-emp0)/2,E[i].type,E[i].u,E[i].v,E[i].sur,E[i+1].sur);
        }
    }
    fclose(out);
}