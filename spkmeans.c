#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
//#include "spkmeans.h"
//#include "kmeans.c"



int getlines(char arr[]) {  
    
    char ch;
    int lines =0;
    int commas_num=0;
    

    FILE *txt = NULL;
    txt = fopen(arr,"r");
    if(txt == NULL){
        printf("Invalid Input! get");
        exit(1);
    } 
    
    while (!feof(txt)){
        ch = fgetc(txt);
        if(lines ==0 && ch ==','){
            commas_num++;
        }
        if (ch == '\n'){
            lines++;
        }
    }
    fclose(txt);
    return lines;
}

/**
 * @brief Get the Dimention of the inserted vectors 
 * 
 * @param arr tne input file 
 * @return int - the vectors dimention
 */

int getDimention(char arr[]){ 
    FILE *txt = NULL;
    char ch;
    int lines =0;
    int commas_num=0;
    txt = fopen(arr,"r");

    if(txt == NULL){
        printf("Invalid Input! size");
        exit(1);
    } 

    while (!feof(txt)){
        ch = fgetc(txt);
        if(lines ==0 && ch ==','){
            commas_num++;
        }
        if (ch == '\n'){
            lines++;
        }
    }
    fclose(txt);
    return commas_num+1;
}

/**
 * @brief read the input file and create the matrix represents the n vectors 
 * 
 * @param input file directory 
 * @param n number of vectors 
 * @param d the vector's dimention
 * @return double** - the matrix represents the n vectors 
 */
 
double** formInputToMatrix(char* input, int n, int d){ //Guy
    int i = 0;
    int j = 0;
    int cnt = 0;
    FILE *txt = NULL;
    double **vectors = (double**) calloc (n,sizeof(double*));
    double *array =  calloc (n*d,sizeof(double));
    txt = fopen(input,"r");


    if (vectors == NULL || txt == NULL || array == NULL){
         printf("An Error Has Occurred");
         exit(1);
    }
    
    while (fscanf(txt, "%lf,", &array[i++])!=EOF){
        
    }
    fclose(txt);
    for(i = 0; i < n; i++){
        vectors[i]=(double*) calloc(d,sizeof(double));
        if (vectors[i] == NULL){
         printf("An Error Has Occurred");
         exit(1);
        }
    }
    for (i = 0; i < n; i++){
        for(j = 0; j < d; j++){
            vectors[i][j] = array[cnt];
            cnt ++;
        }
    } 
    return vectors;
}

void printMatrix(double** matrix, int n, int d){

    for(int i = 0; i < n; i++){
        for (int j =0; j < d; j++){
            printf("%f",matrix[i][j]);
            printf("%s"," ");
        }
        printf("%s\n","");
    }
}

double calcWeight(double* a, double* b, int d)
{
    double sum = 0, norm=0;
    int i=0;
    for (i=0;i<d;i++)
    {
        sum += (a[i]-b[i])*(a[i]-b[i]);
    }
    norm = pow(sum,0.5);
    norm = -norm/2;
    return exp(norm);
}

double** formWeightedMatrix(double** vectorsMatrix, int n, int d)// Yair

{
    double** weightedMatrix =(double **) calloc (n,sizeof(double*));
    assert(weightedMatrix!=NULL && "An Error Has Occured");
    double weight;
    int i=0, row=0, col=0;
    for (i=0; i<n; i++){
        weightedMatrix[i] = calloc (n,sizeof(double));
        assert( weightedMatrix[i]!=NULL && "An Error Has Occured");
    }
    for(row=0;row<n;row++)
    {
        for(col=row+1;col<n;col++)
        {
            weight = calcWeight(vectorsMatrix[row],vectorsMatrix[col],d);
            weightedMatrix[row][col] = weight;
            weightedMatrix[col][row] = weight;
        }
    }
    return weightedMatrix;


}
/**
 * @brief This function generates a degree matrix by all the '1' numbers of each row.
 * There are two cases: Regularmode (in )
 * 
 * 
 * @param weightedMatrix 
 * @param n - number of vertices
 * @param isRegularMode - If isRegularMode ==1 , it means that we want to compute the sum of each row, and if isRegularMode ==0 , it means 
 * that we in sqrt mode.
 * @return double** 
 */

double** formDegreeMatrix (double** weightedMatrix, int n, int isRegularMode){ // Yair
    double * degreeMatrix =(double *) calloc (n,sizeof(double*));
    assert(degreeMatrix!=NULL && "An Error Has Occured");
    int k=0, row,col;
    double sum;
    for (k=0; k<n; k++){
        degreeMatrix[k] = calloc(n,sizeof(double));
        assert( degreeMatrix[k]!=NULL && "An Error Has Occured");
    }

    for (row=0;row<n;row++){
        sum = 0;
        for (col=0;col<n;col++){
            sum += weightedMatrix[row][col];
        }
        if (isRegularMode ==1){
             degreeMatrix[row][row] =(sum);
        }
        else {
            degreeMatrix[row][row] = 1/sqrt(sum);
        }
       
    }
    return degreeMatrix;

}

double** multiplyMatrix(double** aMatrix, double** bMatrix, int n){// Yair

}

double** formLnormMatrix (double** weightedMatrix, int n){ // Yair

}

void formRotaionMatrix(double** P ,double** lNormMatrix, int n){
    // find the largest element off the diaginal
    int i;
    int j;
    int maxRow = 0;
    int maxCol = 0;
    int maxValue;
    double s;
    double c; 
    double t;
    double theta;
    double sign;

    for (i = 0; i < n; i++){ // change for n / 2;
        for(j = 0; j < n; j++){
            if((fabs(lNormMatrix[i][j]) > fabs(lNormMatrix[maxRow][maxCol])) && i != j ){
                maxRow = i;
                maxCol = j;
            }
        }
    }

    theta = (lNormMatrix[maxCol][maxCol] - lNormMatrix[maxRow][maxRow]) / (2 * lNormMatrix[maxRow][maxCol]);
    printf("%f\n",theta);
    sign = theta < 0 ? : 1;
    t = sign / (fabs(theta) + sqrt(pow(theta,2) + 1));
    c = 1 / (sqrt(pow(t, 2) + 1));
    s = t * c;

    for (i = 0; i < n; i++){ // fill in the new  rotation matrix 
        for(j = 0; j < n; j++){
            if ((i == maxRow && j == maxRow) || (i == maxCol && j == maxCol)) // should be complete
                P[i][j] = c;
            else if (i == j)
                P[i][j] = 1;
            else if (i == maxRow && j == maxCol) // seperate foe c 
                P[i][j] = s;
            else if (i == maxCol && j == maxRow)
                P[i][j] = (-1) * s;
            else
                P[i][j] = 0;

        }
    }

    //printMatrix(P,n,n); 
}

int isDiagonal(double** A, int n){
    int i = 0;
    int j = 0;

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(j != i && A[i][j] != 0){
                return 0;
            }

        }
    }
    return 1;

}

void formIdentityMatrix(double** V, int n){ //Guy

    int i;
    for( i = 0; i < n; i++){
        V[i][i] = 1;
    }

}

double ** getTransposeMatrix(double** matrix, int n){
    int i;
    int j;
    double temp;

    for(i = 0; i < n; i++){
        for(j = i + 1; j < n; j++){
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
    return matrix;
}

double** getSimilarMatrix(double** A, double** P,int n){

    double** Pt;
    A = multiplyMatrix(A,P,n);
    P = getTransposeMatrix(P,n);
    return (multiplyMatrix(Pt,A,n));

}

double** jaccobiAlgorithm (double** lNormMatrix , int n){ // Guy - get All the eigenValues

    int i;
    int j;
    double** A = lNormMatrix;

    double** eigenVectors = (double **) calloc (n,sizeof(double*));
    assert(eigenVectors!=NULL && "An Error Has Occured");

    double** P = (double **) calloc (n,sizeof(double*));
    assert(P!=NULL && "An Error Has Occured");

    double** V = (double **) calloc (n,sizeof(double*));
    assert(V!=NULL && "An Error Has Occured");

    for (i=0; i<n; i++){
        P[i] = calloc (n,sizeof(double));
        assert(P[i]!=NULL && "An Error Has Occured");
        eigenVectors[i] = calloc (n,sizeof(double));
        assert(eigenVectors[i]!=NULL && "An Error Has Occured");
        V[i] = calloc (n,sizeof(double));
        assert(P[i]!=NULL && "An Error Has Occured");
    }

    formIdentityMatrix(V,n);
    return V;

    int iteration = 100;

    while(isDiagonal(A, n) == 0 || iteration == 0){ // while A is not Diagonal
        formRotaionMatrix(P,A,n);
        V = multiplyMatrix(V,P,n);
        A = getSimilarMatrix(A,P,n);
        printMatrix(A,n,n);

        iteration --;
    }

    // find the eigenVectors

    return V;
    

}

double** getUmatrix(double** Kvectors){ // Guy 

}



int main(int argc, char *argv[]){ // Guy

    int k = atoi(argv[1]);
    char* input ;
    char* flow;
    int d;
    int n;
    int sizeofvector;
    double ** observationsMatrix;
    
    flow = argv[2];
    input  = argv[3];
    
    n = getlines(input);
    d = getDimention(input);
    observationsMatrix = formInputToMatrix(input, n,d);
    return 0;

}




