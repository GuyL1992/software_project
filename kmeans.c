#define PY_SSIZE_T_CLEAN
// #include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include<math.h>


// static PyObject* k_means_api(PyObject *self, PyObject *args);
// static void fill_matrix_of_points(PyObject* matrix_py, double** matrix_c, int n, int d);
double calcdelta (double *a, double *b, int d);
void addtwovector (double *a, double *b, int d);
void createCentroids(int k, int maxiter, double **init, double **p, int n, int d, double eps);
// static PyObject* kmeans_c(int k, int max_iter, int n, int d, PyObject* vector_py, int epsilon, PyObject* init_centroids_py);



// static PyObject* k_means_api(PyObject *self, PyObject *args)
// {
//     PyObject* vectors_py;
//     PyObject* init_centroids_py;
//     int n;
//     int d;
//     int k;
//     double epsilon;
//     int max_iter;
    
//     if(!PyArg_ParseTuple(args,"iiiidOO",&n,&d,&k,&max_iter,&epsilon,&vectors_py,&init_centroids_py)){ 
//         return NULL;
//     }

//     return Py_BuildValue("O", kmeans_c(k, max_iter,n,d,vectors_py,epsilon,init_centroids_py));
// }

// PyObject* kmeans_c(int k, int max_iter, int n, int d, PyObject* vectors_py, int epsilon, PyObject* init_centroids_py){
//     double** vectors; 
//     double** init_centroids; 
//     int i;
//     int j;
//     PyObject* k_centroids; 
//     PyObject* curr_centroid;

//     vectors = (double **)calloc(n, sizeof(double*));
//     if(vectors == NULL){
//         printf("An Error Has Occurred");
//         exit(1);
//         }

//     for(i = 0; i < n; i++){
//         vectors[i] = calloc(d,sizeof(double));
//         if(vectors[i] == NULL){
//         printf("An Error Has Occurred");
//         exit(1);
//         }
//     }

//     init_centroids = (double **)calloc(k,sizeof(double*));
//     if(init_centroids == NULL){
//         printf("An Error Has Occurred");
//         exit(1);
//         }

//     for(i = 0; i < k; i++){
//         init_centroids[i] = calloc(d,sizeof(double));
//         if(init_centroids[i] == NULL){
//         printf("An Error Has Occurred");
//         exit(1);
//         }
//     }

//     fill_matrix_of_points(vectors_py,vectors,n,d);
//     fill_matrix_of_points(init_centroids_py,init_centroids,k,d);

//     createCentroids(k,max_iter,init_centroids,vectors,n,d,epsilon);
//     k_centroids = PyList_New(k);

//     for(i = 0; i < k; i++){
//         curr_centroid = PyList_New(d);
//         for(j = 0; j < d; j++){
//             PyList_SetItem(curr_centroid,j,PyFloat_FromDouble(init_centroids[i][j]));
//         }
//         PyList_SetItem(k_centroids,i,curr_centroid);
//         free(init_centroids[i]);
//     }
//     free(init_centroids);

//     for(i = 0; i < n; i++){
//        free(vectors[i]);
//     }
//     free(vectors);

//     return k_centroids;
    
// } 



// static void fill_matrix_of_points(PyObject* matrix_py, double** matrix_c, int n, int d){
//     int i;
//     int j;
//     PyObject* vector_py;
//     for (i = 0; i < n; i++){
//         vector_py = PyList_GET_ITEM(matrix_py,i);
//         if(vector_py == NULL){
//         printf("An Error Has Occurred");
//         exit(1);
//         }
//         for(j = 0; j < d; j++){
//             matrix_c[i][j] = PyFloat_AsDouble(PyList_GET_ITEM(vector_py,j));
//         }        
//     }
// }

double calcdelta (double *a, double *b, int d){ 

    double sum = 0.0;
    int i = 0;

    for(i = 0 ; i < d; ++i)
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    
    return pow(sum,0.5);  
}
   
void addtwovector (double *a, double *b, int d){  
    int i = 0;
    for(i = 0 ; i < d; ++i){
        a[i] += b[i];
    }
}

void calculateKCentroids(int k, int maxiter,double **init, double** observations,int n, int d, double eps){ 

    const double epsilon = eps;
    int norm = 1; 
    int cnt = 0;
    int i = 0;
    int j = 0;
    int t = 0;

    double** addvectors = (double **)calloc(k, sizeof(double*));
    int* sizes = (int*)calloc(k, sizeof(int));

    printf("%s", "calculating cectroids");

    if (addvectors == NULL || sizes == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    
    for (i = 0; i < k; i++){
        addvectors[i] = (double *)calloc(d, sizeof(double));
        if (addvectors[i] == NULL){
         printf("An Error Has Occurred");
         exit(1);
        }
        for (j = 0; j < d; j++){
            addvectors[i][j] = 0;
        }
    }

    while(maxiter > 0 && norm == 1){
        
        for (i = 0; i < n; i++){
            double min_delta = 1000000.0;
            int c_index = -1;

            for (j = 0; j < k ; j++){
                double temp = calcdelta(observations[i],init[j], d);
                if (temp <= min_delta){
                    c_index = j;
                    min_delta = temp;
                }   
            }  
            addtwovector(addvectors[c_index],observations[i],d);
            sizes[c_index]++;
        }

        cnt = 0;
        for (i = 0; i < k; i++){
            double *prev = (double *)calloc(d, sizeof(double));
            double dis;
            if (prev == NULL){
                printf("An Error Has Occurred");
                exit(1);
            }

            for (j = 0; j < d; j++){
                prev[j] = init[i][j];
            }

            for (t = 0; t < d; t++){
                init[i][t] = addvectors[i][t] / sizes[i];
            }
            
            dis = calcdelta(prev,init[i],d);
            free(prev); 

            if (dis < epsilon)
                cnt += 1;
        }

        if (cnt == k)
            norm = 0;

        for (i = 0; i< k; i++ ) {
            sizes[i] = 0;
        }

         for (i = 0; i < k; i++){
             for (j = 0; j < d; j++)
             addvectors[i][j]=0;
         }

        maxiter --;
        
  
    }
    free(sizes);

    for (i = 0; i < k; i++){
        free(addvectors[i]);
    }
    free(addvectors);
}

// static PyMethodDef _capiMethods[] = {
//     {"fit", (PyCFunction) k_means_api, METH_VARARGS, PyDoc_STR("function to calculate and return k Clusters")},
//     {NULL,NULL,0,NULL}
// };

// static struct PyModuleDef _moduledef = {
//     PyModuleDef_HEAD_INIT,
//     "mykmeanssp",
//     NULL,
//     -1,
//     _capiMethods
// };

// PyMODINIT_FUNC
// PyInit_mykmeanssp(void)
// {
//     PyObject *m;
//     m = PyModule_Create(&_moduledef);
//     if(!m){
//         return NULL;
//     }
//     return m;
// }




