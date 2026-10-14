#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define frand (rand()/((double) RAND_MAX+1))
#define N 2000
#define k 10
#define Npasos 1000
#define k_4 40
//#define DEBUG



void ErdosReyni (int A[][k_4]){
    double p;
    int i = 0, j = 0,l;
    double r = 0;
    p=(double)k/(N-1);
    for (i = 0; i<N; i++){
        l=0;
        for (j = 0; j<i; j++){
            r = frand;
            if (r<p){
                while (A[i][l]!=0)
                    l++;
                A[i][l]=j+1;
                l=0;
                while (A[j][l]!=0)
                    l++;
                A[j][l]=i+1;
                l=0;
            }
        }
    }
}
double Conexa (int B[][k_4]){
    int i,j,n[N],m[N],aux,l;
    aux=1;
    double nnodos;
    nnodos=N;
    for(i=0;i<N;i++){ //Inicializamos el array
        n[i]=0;
        m[i]=0;
    }
    n[0]=1;
    while(aux==1){
        for(i=0;i<N;i++)
            for(j=0;B[i][j]!=0;j++)
                n[B[i][j]-1]=1;
        for(i=0,aux=0;i<N;i++){
            if(n[i]!=m[i]){
                aux=1;
                m[i]=n[i];
            }
        }
    }
    for(i=0;i<N;i++){
        if(n[i]==0){
            nnodos--;
            for(l=i;l<N-1;l++){
                j=0;

                while (j<k_4){
                    B[l][j]=B[l+1][j];
                    B[l+1][j]=0;
                    j++;
                }
            }
        }
    }
    return nnodos;
}
void WattsStrogatz (double p, int A[][k_4]){
    int i,j;
    int e;
    double r;
    e=0;
    for (i=0; i<N; i++){
        for (j=0;j<k_4; j++){
            A[i][j]=0;
        }
    }
    for (i=0; i<N;i++){
        e=0;
        for (j=i-k/2;j<i;j++){
            A[i][e] = j+1;
            e++;
        }
        for(j=i+1;j<=k/2+i;j++){
            A[i][e] = j+1;
            e++;
        }

    }
    for (i=0;i<N;i++){
        for (j=0;j<k;j++){
            if (A[i][j]<=0){
                A[i][j]+=N;
            }
            if (A[i][j]>N){
                A[i][j]-=N;
            }
            //printf("%d \t" , A[i][j]);
        }
        //printf("\n");

    }
    for (i=0; i<N;i++){
        for (j=0;A[i][j]!=0;j++){
            r=frand;
            if (p>r){
                A[i][j]=(int)N*frand+1;
                if (A[i][j]==i+1)
                    i--;
            }
               //printf("%d \t" , A[i][j]);
        }
        //printf("\n");
    }
}

void f2(int B[][k_4],double *x, double *derivada, double *incremento, double div,double beta, double delta, double K){
    int i,j,l, auxiliar;
    double w[N][k];
    double xi,xj;
    double suma,suma_2;
    /*for (i=0;i<N;i++){
     for (j=0;j<k;j++){
     w[i][j]=-1000.0;
     }
     }*/
    for (i=0;i<N;i++){
        xi=x[i]+incremento[i]/div;
        for(j=0;j<k;j++){
            xj=x[B[i][j]-1]+incremento[j]/div;
            auxiliar=-1;
            if (i>B[i][j]-1){
                for (l=0;l<k;l++){
                    if(B[B[i][j]-1][l]==i+1)
                        auxiliar=l;
                }
                if(auxiliar==-1){
                    w[i][j]=pow((fabs(xi-xj))+delta,-beta);
                }
                if (auxiliar>-1){
                    w[i][j]=w[B[i][j]-1][auxiliar];
                }}
            if (i<=B[i][j]-1){
                w[i][j]=pow((fabs(xi-xj))+delta,-beta);
            }
        }
    }
}
void RungeKutta (FILE *f,int npasos, double dt, double *x,int B[][k_4],double beta, double delta,double K){
    int i,j;
    double cero[N];
    for (i=0;i<N;i++){
        cero[i]=0;
    }
    double J[4][N];
    for (i=0; i<npasos;i++){
        f2(B,x,J[0],cero,1.0/dt,beta, delta, K);
        f2(B,x,J[1],J[0],2.0/dt,beta, delta, K);
        f2(B,x,J[2],J[1],2.0/dt,beta, delta, K);
        f2(B,x,J[3],J[2],1.0/dt,beta, delta, K);
        for(j=0;j<N;j++){
            x[j]+=(J[0][j]+2*J[1][j]+2*J[2][j]+J[3][j])*dt/6;
        }
    }
}
double med(const double *x){
    int i;
    double suma;
    for (i=0,suma=0;i<N;i++)
        suma+=x[i];
    return suma/N;
}
double absmed(const double *x){
    return fabs(med(x));
}
double desvest (const double *x){
    int i;
    double suma,media;
    media=med(x);
    for(i=0,suma=0;i<N;i++)
        suma+=(x[i]-media)*(x[i]-media);
    return sqrt(suma/N);

}
int main (void){
    srand(time(NULL));
    int B[N][k_4];
    double x[N];
    int i,j,l;
    int nnodos;
    double deltat,delta;
    double beta,K;
    double fraccion[10][2];
    FILE *f;
    f=fopen("C:\\Users\\ricar\\Desktop\\CaosTr\\Códigos\8.txt","a");
    deltat=1;
    K=15;
    beta=0.8;
    delta=0.002*K;
    nnodos=0;
    for(i=0;i<50;i++){
        for(j=0;j<N;j++){
            x[j]=(-1+2*frand)*K;
            for(l=0;l<k_4;l++)
                B[j][l]=0;
        }
        ErdosReyni(B);
        RungeKutta(f, Npasos, deltat, x, B, beta, delta, K);
        fraccion[i][0]=absmed(x);
        fraccion[i][1]=desvest(x);
        //printf("%lf %lf\n",absmed(x),desvest(x));
        fprintf(f, "%d\t%lf\t%lf\n",i,fraccion[i][0],fraccion[i][1]);
        printf("%d\n",i);
    }
    return 0;
}
