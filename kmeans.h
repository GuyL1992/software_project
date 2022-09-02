static PyObject* k_means_api(PyObject *self, PyObject *args);
static void fill_matrix_of_points(PyObject* matrix_py, double** matrix_c, int n, int d);
double calcdelta (double *a, double *b, int d);
void addtwovector (double *a, double *b, int d);
void calculateKCentroids(int k, int maxiter, double **init, double **p, int n, int d, double eps);
static PyObject* kmeans_c(int k, int max_iter, int n, int d, PyObject* vector_py, int epsilon, PyObject* init_centroids_py);
