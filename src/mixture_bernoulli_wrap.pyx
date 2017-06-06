from libc.stdlib cimport malloc, free

cdef extern from "mixture_bernoulli.h":
    int mixtureBernoull(unsigned int *sample_X, int n_samples, int n_dimentions, int n_components, double **means)
    
def mixtureBernoullWrap(sample_X, n_components):
    cdef int n_samples = sample_X.shape[0]
    cdef int n_dimentions = sample_X.shape[1]
    cdef double **means = <double**> malloc(sizeof(double*) * n_components)
    cdef unsigned int *c_sample_X = <unsigned int*> malloc(sizeof(unsigned int) * n_samples * n_dimentions)
    cdef int i, j
    
    for i in range(n_components):
        means[i] = <double*> malloc(sizeof(double) * n_dimentions)

    for i in range(n_samples):
        for j in range(n_dimentions):
            c_sample_X[i * n_dimentions + j] = sample_X[i, j]
            
    mixtureBernoull(c_sample_X, n_samples, n_dimentions, n_components, means)
    
    result = []
    for i in range(n_components):
        row = []
        for j in range(n_dimentions):
            row.append(means[i][j])
        result.append(row)
        
    return result
