#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
//#include "spkmeans.h"
//#include "kmeans.c"

int convertStringIntoGoalEnum(char* UserGoal)
{
    if (!strcmp(UserGoal,"spk"))
    {
        return 0;
    } else if (!strcmp(UserGoal,"wam")){
        return 1;
    } else if (!strcmp(UserGoal, "ddg")){
        return 2;
    } else if (!strcmp(UserGoal,"lnorm")){
        return 3;
    } else if (!strcmp(UserGoal,"jacobi")){
        return 4;
    } else 
    {
        return 5;
    };
}
double** allocationMatrix(int n, int d);
void freeMatrix(double** matrix,int n);

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


int checkNegativeZero(double value)
{
    if (value>-0.00005 && value<0){ return 1;}
    return 0;
}

void printMatrix2(double** matrix, int a, int b) {
    int i,j;
    for (i = 0; i < a; i++)
    {
        for (j = 0; j <  b - 1; j++)
        {
            if (checkNegativeZero(matrix[i][j])){
                printf("0.0000,");
            } else{
                printf("%.4f,", matrix[i][j]);
            }
            
        }
        if (checkNegativeZero(matrix[i][b-1])){
            printf("0.0000");
        } else{
            printf("%.4f", matrix[i][b-1]);
        }
        putchar('\n');
    }
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
    double sum = 0, norm=0, currDiff;
    int i=0;
    for (i=0;i<d;i++)
    {
        currDiff = (a[i]-b[i]) * (a[i]-b[i]);
        sum = sum + currDiff;
    }
    norm = pow(sum,0.5);
    norm = -norm/2;
    return exp(norm);
}

void formWeightedMatrix(double** weightedMatrix,double** vectorsMatrix, int n, int d)// Yairr

{
    double weight;
    int row;
    int col;
    
    for(row=0;row<n;row++)
    {
        for(col=row+1;col<n;col++)
        {
            weight = calcWeight(vectorsMatrix[row],vectorsMatrix[col],d);
            weightedMatrix[row][col] = weight;
            weightedMatrix[col][row] = weight;
        }
    }
    return;


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

void formDegreeMatrix (double** degreeMatrix, double** weightedMatrix, int n, int isRegularMode){ // Yair
    int k=0, row,col;
    double sum;

    for (row=0;row<n;row++){
        sum = 0;
        for (col=0;col<n;col++){
            sum += weightedMatrix[row][col];
        }
        if (isRegularMode == 1){
             degreeMatrix[row][row] =(sum);
        }
        else {
            degreeMatrix[row][row] = 1/sqrt(sum);
        }
       
    }
    return;

}

void multiplyMatrix(double** result, double** aMatrix, double** bMatrix, int n){// Yair
    
    int k,i,j;
    double sum;

    for( i = 0; i < n; i++){
        for(j = 0; j < n; j ++){
            result[i][j] = 0;
        }
    }

    
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            sum = 0;
            for(k=0;k<n;k++)
            {
                result[i][j]+=aMatrix[i][k] * bMatrix[k][j];
            }
        }
    }
    return;
}

void formLnormMatrix (double** lNormMatrix, double** weightedMatrix,double** degreeSqrtMatrix, int n){ // Yair

    double ** temp = allocationMatrix(n,n);
    multiplyMatrix(temp, degreeSqrtMatrix, weightedMatrix,n);
    multiplyMatrix(lNormMatrix, temp, degreeSqrtMatrix,n);

int i=0,j=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (i!=j)
            {
                lNormMatrix[i][j] = -lNormMatrix[i][j];
                
            }
            else
            {
                lNormMatrix[i][i] = 1 - lNormMatrix[i][i];
            }
        }
    }
    freeMatrix(temp,n);
    return;
}

void formRotaionMatrix(double** P ,double** A,double** matrixV, int n){
    // find the largest element off the diaginal
    int i;
    int j;
    int r;
    int maxRow;
    int maxCol;
    double maxValue = 0;
    double s;
    double c; 
    double t;
    double theta;
    double sign;



    for (i = 0; i < n; i++){ // change for n / 2;
        for(j = 0; j < n; j++){
            if(fabs(A[i][j]) > fabs(maxValue) && i != j ){
                maxRow = i;
                maxCol = j;
                maxValue = A[maxRow][maxCol];
            }
        }
    }

    // theta = (A[maxCol][maxCol] - A[maxRow][maxRow]) / (2 * A[maxRow][maxCol]);
    // sign = theta < 0 ? : 1;
    // t = sign / (fabs(theta) + sqrt(pow(theta,2) + 1));
    // c = 1 / (sqrt(pow(t, 2) + 1));
    // s = t * c;

    double tetha, signTetha, absTetha;
    
    tetha = (A[j][j]-A[i][i])/(2*A[i][j]);
    (tetha<0) ? (signTetha = -1) : (signTetha=1);
    (tetha<0) ? (absTetha = -tetha) : (absTetha = tetha);
    t = signTetha / (absTetha + pow((pow(tetha,2)+1),0.5));
    c = 1/(pow((pow(t,2)+1),0.5));
    s = t*(c);

   

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

    double *Icol, *Jcol;
    i = maxRow;
    j = maxCol;
    double** matrixA = A;

    Icol = calloc(n,sizeof(double));
    assert(Icol!=NULL && "An Error Has Occured");
    Jcol = calloc(n,sizeof(double));
    assert(Jcol!=NULL && "An Error Has Occured");

    for (r=0;r<n;r++)
    {
        Icol[r] = matrixA[r][i];
        Jcol[r] = matrixA[r][j];
    }
    for(r=0; r<n;r++)
    {
        if(r!=i && r!=j)
        {
            matrixA[r][i] = c*Icol[r] - s*Jcol[r];
            matrixA[i][r] = matrixA[r][i];
            matrixA[r][j] = c*Jcol[r] + s* Icol[r];
            matrixA[j][r] = matrixA[r][j];
        }
    }

    matrixA[i][i] = pow(c,2)*Icol[i] + pow(s,2)*Jcol[j] - 2 * s * c * Jcol[i];
    matrixA[j][j] = pow(s,2) *Icol[i] + pow(c,2) * Jcol[j] + 2 * s * c * Jcol[i];
    matrixA[i][j] = 0;
    matrixA[j][i] = 0;
    
    free(Icol);
    free(Jcol);

    Icol = calloc(n,sizeof(double));
    assert(Icol!=NULL && "An Error Has Occured");
    Jcol = calloc(n,sizeof(double));
    assert(Jcol!=NULL && "An Error Has Occured");

    for (r=0;r<n;r++)
    {
        Icol[r] = matrixV[r][i];
        Jcol[r] = matrixV[r][j];
    }
    for (r=0;r<n;r++)
    {
        matrixV[r][i] = c*Icol[r] - s*Jcol[r];
        matrixV[r][j] = s*Icol[r] + c*Jcol[r];
    }
    free(Icol);
    free(Jcol);

    return;
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

void getSimilarMatrix(double** A, double** P,double** temp, int n){

    double** Pt;
    multiplyMatrix(temp,A,P,n);
    Pt = getTransposeMatrix(P,n);
    multiplyMatrix(A,Pt,temp,n);
    
    return ;

}

double getOff(double ** A, int n){
    double offA = 0;
    int i = 0;
    int j = 0;

    for(i=0;i<n;i++)
    {
        for(j=i+1;j<n;j++)
        {
            offA += 2 * pow(A[i][j],2);
        }
    }

    return offA;
}

void copyMatrix (double** A, double** B, int n){
    int i;
    int j;

    for(i = 0; i < n; i++){
        for (j = 0; j < n; j ++){
            A[i][j] = B[i][j];
        }
    }

    return;
}




void jaccobiAlgorithm (double** eigenVectors,double** observationsMatrix , int n){ // Guy - get All the eigenValues

    int i;
    int j;
    int iteration = 0;
    double** A = observationsMatrix;
    double** V = eigenVectors;
    double offA;
    double offAtag = getOff(A,n);
    double epsilon = pow(10,-15);
    int stopCondition = 1;
    double** P = allocationMatrix(n,n);
    double** helperPointer;
    double** temp = allocationMatrix(n,n);
    double** Atag = allocationMatrix(n,n);



    formIdentityMatrix(V,n); // first: V equal to the Identity matrix 

    do{
        offA = offAtag;
        formRotaionMatrix(P,A,V,n);
        // multiplyMatrix(temp,V,P,n);
        // copyMatrix(V,temp,n);
        // getSimilarMatrix(A,P,temp,n);
        offAtag = getOff(A,n);
        iteration ++;

       if ((offA-offAtag)<epsilon)
        {
            stopCondition = 0;
        }
        if (offAtag==0) {break;}
        



    } while(stopCondition && iteration<100);

    printMatrix2(V,n,n);
    
    freeMatrix(temp,n);
    freeMatrix(P,n);

    return;
 
}

void getUmatrix(double** Kvectors){ // Guy 
    return;

}

double** allocationMatrix(int n, int d){

    double** allocatedMatrix =(double **) calloc (n,sizeof(double*));
    assert(allocatedMatrix!=NULL && "An Error Has Occured");
    int i=0, row=0, col=0;

    for (i=0; i<n; i++){
        allocatedMatrix[i] = calloc (n,sizeof(double));
        assert(allocatedMatrix[i]!=NULL && "An Error Has Occured");
    }

    return allocatedMatrix;
}

void freeMatrix(double** matrix,int n){

    int i=0;

    for (i=0; i<n; i++){
        free(matrix[i]); 
    }

    free(matrix);

    return;
}
   



void process_wam(double** observationMatrix, int n, int d){
    double** weightedMatrix = allocationMatrix(n,n);
    formWeightedMatrix(weightedMatrix,observationMatrix,n,d);
    printMatrix2(weightedMatrix,n,n);
    freeMatrix(weightedMatrix,n);
}

void process_ddg (double** observationMatrix, int n, int d){
    double** weightedMatrix = allocationMatrix(n,n);
    double** degreeMatrix = allocationMatrix(n,n);
    int regularMode = 1;
    formWeightedMatrix(weightedMatrix,observationMatrix,n,d);
    formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    freeMatrix(weightedMatrix,n);
    printMatrix2(degreeMatrix,n,n);
    freeMatrix(degreeMatrix,n);

}

void process_lnorm (double** observationsMatrix, int n, int d){
    int regularMode = 0;
    double** weightedMatrix = allocationMatrix(n,n);
    double** degreeMatrix = allocationMatrix(n,n);
    double** lNormMatrix = allocationMatrix(n,n);
    formWeightedMatrix(weightedMatrix,observationsMatrix,n,d);
    formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    formLnormMatrix(lNormMatrix, weightedMatrix,degreeMatrix,n);

    freeMatrix(weightedMatrix,n);
    freeMatrix(degreeMatrix,n);

    printMatrix2(lNormMatrix,n,n);
    freeMatrix(lNormMatrix,n);

}

void jaccobi_proccess(double** observationsMatrix, int n, int d){

    int regularMode = 0;

    // double** weightedMatrix = allocationMatrix(n,n);
    // double** degreeMatrix = allocationMatrix(n,n);
    // double** lNormMatrix = allocationMatrix(n,n);
    // formWeightedMatrix(weightedMatrix,observation_matrix,n,d);
    // formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    // formLnormMatrix(lNormMatrix, weightedMatrix,degreeMatrix,n);
    // freeMatrix(weightedMatrix,n);
    // freeMatrix(degreeMatrix,n);

    double** eigenVectors = allocationMatrix(n,n);
    // double* eigenValues = calloc(n, sizeof(double));

    jaccobiAlgorithm(eigenVectors,observationsMatrix,n);

    printMatrix2(eigenVectors,n,n);
    // free(eigenValues);
    freeMatrix(eigenVectors,n);

    return;

}

     


/**
 * @brief This function determine the number of cluskers k. k is the max gap between two 
 * following eingevalues, until half of the values.
 * 
 * @param eigenValues 
 * @param len len of eingevalues.
 * @return int - k 
 */

 int cmpfunc (const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
}

int TheEigengapHeuristic(double* eigenValues, int lenOfArr) {

    int index=0;
    double maxDelta = 0;
    double delta = 0;
    int i;
    
    qsort(eigenValues,lenOfArr, sizeof(double), cmpfunc);
    for(i=1; i<=(lenOfArr/2);i++)
    {
        delta = eigenValues[i]-eigenValues[i-1];
        if (delta > maxDelta){
            maxDelta = delta;
            index = i;
        }
    }
    
    return index;
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
    //jaccobi_proccess(observationsMatrix,n,d);

    switch(convertStringIntoGoalEnum(flow)){

        case 0:
            break;

        case 1:
            process_wam(observationsMatrix,n,d);
            break;

        case 2:
            process_ddg(observationsMatrix,n,d);
            break;

        case 3:
            process_lnorm(observationsMatrix,n,d);
            break;

        case 4:
            jaccobi_proccess(observationsMatrix,n,d);
            break;
        
        default:
            break;


    }//

    freeMatrix(observationsMatrix,n);
    
    
    return 0;

}




