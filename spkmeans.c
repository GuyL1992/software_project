#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
//#include "spkmeans.h"
//#include "kmeans.c"


/
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

double calculateWeight(double* a, double* b, int d)
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
    double * weightedMatrix =(double *) calloc (n,sizeof(double*));
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
        degreeMatrix[k] = calloc (n,sizeof(double));
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

double** formLnormMatrix (double** weightedMatrix, int n){ // Guy

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
    printMatrix(observationsMatrix,n,d);

    return 0;

}




