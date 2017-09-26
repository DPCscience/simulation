//https://www.daniweb.com/programming/software-development/threads/320739/matrix-struct-and-adding-matrices
/* Allocating mem for m.mat */
mat_temp = (float**) malloc(a.rows * sizeof(float*));
for (i = 0; i < a.rows; i++){
    mat_temp[i] = (float*) malloc(a.cols * sizeof(float));
}
m.rows = a.rows;
m.cols = m.cols;
m.mat = mat_temp;


/* Structures */
typedef struct matrix {
    int rows;
    int cols;
    float **mat;
} Matrix;


Matrix matrix_add(Matrix a, Matrix b){
    Matrix m;
    int i, j;
    for (i=0; i<a.rows; i++) {
        for (j=0; j<a.cols; j++) {
            m.mat[i][j] = (a.mat[i][j]+b.mat[i][j]);
        }
    }
    m.rows = a.rows;
    m.cols = a.cols;
    return m;
}