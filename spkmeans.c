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

double** formWeightedMatrix(double** vectorsMatrix, int n)// Yair
{
}

double** formDegreeMatrix (double** weightedMatrix, int n){ // Yair

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




