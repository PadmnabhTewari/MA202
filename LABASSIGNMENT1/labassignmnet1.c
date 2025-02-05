#include<stdio.h>
#include<math.h>
void transpose(double A[3][3],double B[3][3]) 
{ 
    int i, j; 
    for (i = 0; i < 3; i++) 
        for (j = 0; j < 3; j++)
            B[i][j] = A[j][i]; 
} 
void copy(double A[3][3],double B[3][3]){
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            B[i][j]=A[i][j];
        }
    }
}
void jacobii(double J[3][3],double x[3][3],int n)
{
    int a=0,b=0; 
    n=n-1;
    double max=0;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            if(i==j){
                continue;
            }
            else{
                if(J[i][j]>max){
                    max=J[i][j];
                    a=i;
                    b=j;
                }
            }
        }
    }
    double theta=atan((2*J[a][b])/(J[a][a]-J[b][b]));
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            if((i==a && j==b) || (i==b && j==a)){
                x[i][j]=cos(theta);
                continue;
            }
            if(i==a && j==a){
                x[a][a]=sin(theta);
                continue;
            }
            else if(i==b && j==b){
                x[b][b]=-sin(theta);
                continue;
            }
            else if(i==j){
                x[i][i]=1;
            }
            else{
                x[i][j]=0;
            }
        }
    }
    double c[3][3];
    double J2[3][3];
    copy(J,J2);
    transpose(J,J2);
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        {    
            c[i][j]=0;    
            for(int k=0;k<3;k++)    
            {    
                c[i][j]+=J2[i][k]*x[k][j];    
            }    
        }    
    }
    double d[3][3];
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        {    
            d[i][j]=0;    
            for(int k=0;k<3;k++)    
            {    
                d[i][j]+=c[i][k]*J[k][j];    
            }    
        }    
    }
    int c1=0;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            if(i==j){
                continue;
            }
            else{
                if(d[i][j]!=0){
                    c1++;
                }
            }
        }
    }
    if(c1!=0){
        jacobii(d,x,n);
    }
}
int main(){
    double A[3][3];
    printf("enter the matrix");
    double x[3][3];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            scanf("%lf",&A[i][j]);
            x[i][j]=0;
        }
    }
    double B[3][3];
    transpose(A,B);
    double J[3][3];
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        {    
            J[i][j]=0;    
            for(int k=0;k<3;k++)    
            {    
                J[i][j]+=A[i][k]*B[k][j];
            }    
        }    
    }
    int n=10;
    jacobii(J,x,n);
    double sig[3][3];
    printf("The value is %d",n);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            if(i==j){
                sig[i][i]=sqrt(J[i][i]);
            }
            else{
                sig[i][j]=0;
            }
        }
    }
    double J2[3][3];
    double y[3][3];
    transpose(J,J2);
    n=10;
    jacobii(J2,y,n);
    double Vt[3][3];
    transpose(y,Vt);
    double U[3][3];
    copy(x,U);
    double A1[3][3];
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        {    
            A1[i][j]=0;    
            for(int k=0;k<3;k++)    
            {    
                A1[i][j]+=U[i][k]*sig[k][j];
            }    
        }    
    }
    double ans[3][3];
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        {    
            ans[i][j]=0;    
            for(int k=0;k<3;k++)    
            {    
                ans[i][j]+=A1[i][k]*Vt[k][j];
            }    
        }    
    }
    for(int i=0;i<3;i++){    
        for(int j=0;j<3;j++)    
        { 
            printf("%f ",ans[i][j]);
        }
        printf("\n");
    }
    return 0;
}



