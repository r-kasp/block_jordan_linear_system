#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    int n, m, r, s;
    char * filename = nullptr;
    if (!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && 
        sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1)) 
    {
        printf("Usage %s n m r s (file) \n", argv[0]);
        return 1;
    }
    
    if (argc == 6)
        filename = argv[5];
    //Исходная матрица, которую не будем преобразовывать, чтобы сверить ответ
    double *a = new double[n*n];
    if (!a)
    {
        printf("Not enough memory\n");
        return 2;
    }
    //Вектор b
    double *b = new double[n];
    if (!b)
    {
        printf("Not enough memory\n");
        delete [] a;
        return 3;
    }
    //Копия вектора b
    double *c = new double[n];
    if (!c)
    {
        printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        return 4;
    }
    
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] a;
            delete [] b;
            delete [] c;
            return 5;
        }
    }
    else
        init_matrix(a,n,s);
    init_vector(b, a, n);
        
    print_matrix(a,n,n,r);
    
    copy_vector(c, b, n);
    
    double elapsed = clock();
    int res = block_jordan(a, c, n, m);
    double *rsd = new double[n];
    if (!rsd)
    {
    	  printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        delete [] c;
        return 6;
    }
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] a;
            delete [] b;
            delete [] c;
            return 7;
        }
    }
    else
        init_matrix(a,n,s);
    multiply_matrix_and_vector(a, c, rsd, n, n);
    subtract_vectors(rsd, b, n);
    double residual = norma_vector(rsd, n)/norma_vector(b, n);
    elapsed = (clock()-elapsed)/CLOCKS_PER_SEC;
    printf("Result : \n");
    if (res != SUCCESS)
    {
        printf("Can't find solution. The matrix is singular\n");
        printf("Elapsed = %.2f\n", elapsed);
    }
    else
    {
        print_vector(c,n,r);
        printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], residual, elapsed, s, n, m);
    }
    delete [] a;
    delete [] b;
    delete [] rsd;
    delete [] c;
    return 0;
}
