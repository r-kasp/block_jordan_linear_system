#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "matrix.h"

int read_matrix(double *a, int n, char *filename)
{
    FILE *fp;
    int i,j;
    if (!(fp = fopen(filename,"r")))
        return ERROR_OPEN;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (fscanf(fp, "%lf", a+i*n+j) != 1)
            {
                fclose(fp);
                return ERROR_READ;
            }
    fclose(fp);
    return SUCCESS;
}

double f(int s, int n, int i, int j)
{
    i++; j++;
    if (s == 1)
        return n - (i > j ? i : j) + 1;
    if (s == 2)
        return (i > j ? i : j);
    if (s == 3)
        return (i > j ? i-j : j-i);
    if (s == 4)
    {
        return (double)1/(i+j-1); 
    }
    return 1;
}

void init_matrix(double *a, int n, int s)
{
    int i,j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i*n+j] = f(s,n,i,j);
}

void print_matrix(double *a, int m, int n, int r)
{
    double n_max = (n > r ? r : n);
    double m_max = (m > r ? r : m);
    for (int i = 0; i < m_max; i++)
    {
        for (int j = 0; j < n_max; j++)
        {
            printf(" %10.3e", a[i*n+j]);
        }
        printf("\n");
    }
}

double norma_vector(double *b, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
    	double check = fabs(b[i]);
        if (check > res)
        	res = check;
    }
    return res;
}


//horizontal norm
double norma_matrix(double *a, int m, int n)
{
  double res = 0;
  for (int i = 0; i < m; i++)
  {
    double check = 0;
    for (int j = 0; j < n; j++)
    {
      check += fabs(a[i*n + j]);
    }
    if (check > res)
    {
      res = check;
    }
  }
  return res;
}


void multiply_matrix(double *a, double *b, double *c, int m, int n, int k)
{
    for (int i = 0; i < m; i++)
    {
        double * cbuf = c + i * k;
        for (int j = 0; j < k; j++)
          cbuf[j] = 0;
        for (int j = 0; j < n; j++)
        {
          double mult = a[i * n + j];
          double *bbuf = b + j * k;
          for (int q = 0; q < k; q++)
            cbuf[q] += mult * bbuf[q];
        }
    }
}


void multiply_blocks_m(double *a, double *b, double *c, int m)
{
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    int to = (m / 3) * 3;
    //идём до последней строчки которая не поделилась
    for (int i = 0; i < to; i += 3)
    {
      //сначала по обычным блокам
      for (int j = 0; j < to; j += 3)
      {
        s00 = 0; 
        s01 = 0; 
        s02 = 0;
        s10 = 0; 
        s11 = 0; 
        s12 = 0;
        s20 = 0; 
        s21 = 0; 
        s22 = 0;
        for (int q = 0; q < m; q++)
        {
          double a00 = a[i * m + q], a10 = a[(i + 1) * m + q], a20 = a[(i + 2) * m + q],
                 b00 = b[q * m + j], b01 = b[q * m + j + 1], b02 = b[q * m + j + 2];
          s00 += a00 * b00; 
          s01 += a00 * b01; 
          s02 += a00 * b02;
          s10 += a10 * b00; 
          s11 += a10 * b01; 
          s12 += a10 * b02;
          s20 += a20 * b00; 
          s21 += a20 * b01; 
          s22 += a20 * b02;
        }
        c[i * m + j] = s00; 
        c[i * m + j + 1] = s01; 
        c[i * m + j + 2] = s02;
        c[(i + 1) * m + j] = s10; 
        c[(i + 1) * m + j + 1] = s11; 
        c[(i + 1) * m + j + 2] = s12;
        c[(i + 2) * m + j] = s20; 
        c[(i + 2) * m + j + 1] = s21; 
        c[(i + 2) * m + j + 2] = s22;
      }
      //теперь остался один неполный блок
      for (int j = to; j < m; j++)
      {
        s00 = 0; 
        s10 = 0; 
        s20 = 0;
        for (int q = 0; q < m; q++)
        {
          double b00 = b[q * m + j];
          double a00 = a[i * m + q], a10 = a[(i + 1) * m + q], a20 = a[(i + 2) * m + q];
          s00 += a00 * b00;
          s10 += a10 * b00;
          s20 += a20 * b00;
        }
        c[i * m + j] = s00;
        c[(i + 1) * m + j] = s10;
        c[(i + 2) * m + j] = s20;
      }
    } 
    //последняя не поделившаяся строчка
    for (int i = to; i < m; i++)
    {
      //идем до маленького блока
      for (int j = 0; j < to; j += 3)
      {
        s00 = 0; 
        s01 = 0; 
        s02 = 0;
        for (int q = 0; q < m; q++)
        {
          double a00 = a[i * m + q];
          double b00 = b[q * m + j], b01 = b[q * m + j + 1], b02 = b[q * m + j + 2];
          s00 += b00 * a00;
          s01 += b01 * a00;
          s02 += b02 * a00; 
        }
        c[i * m + j] = s00;
        c[i * m + j + 1] = s01;
        c[i * m + j + 2] = s02;
      }
      //маленький блок
      for (int j = to; j < m; j++)
      {
        s00 = 0;
        for (int q = 0; q < m; q++)
          s00 += a[i * m + q] * b[q * m + j];
        c[i * m + j] = s00;
      }
    }
}

/*
void multiply_matrix(double *a, double *b, double *c, int m, int n, int k)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < k; j++)
        {
            double s = 0;
            for (int q = 0; q < n; q++)
                s += a[i*n+q]*b[q*k+j];
            c[i*k+j] = s;
        }
    }
}
*/
void multiply_matrix_and_vector(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
    	  double s = 0;
        for (int j = 0; j < n; j++)
        {
            s += a[i*n + j] * b[j];
        }
        c[i] = s;
    }
}

int solve(int n, double *a, double *x, double *b, double *c)
{
	(void)a;
	(void)x;
	(void)b;
	for (int i = 0; i < n; i++)
		c[i] = (i+1)%2;
	return SUCCESS;
}

void subtract_vectors(double *b, double *c, int n)
{
	for (int i = 0; i < n; i++)
		b[i] -= c[i];
}

void copy_matrix(double *x, double *a, int n)
{
	for (int i = 0; i < n*n; i++)
			x[i] = a[i];
}

void copy_vector(double *c, double *b, int n)
{
	for (int i = 0; i < n; i++)
		c[i] = b[i];
}

void init_vector(double *b, double *a, int n)
{
	for (int i = 0; i < n; i++)
	{
		double s = 0;
		for (int k = 0; k < (n+1)/2; k++)
		{
			s += a[i*n + 2*k];
		}
		b[i] = s;
	}
}

void print_vector(double *b, int n, int r)
{
	double n_max;
  n_max = (n > r ? r : n);
  for (int i = 0; i < n_max; i++)
  	printf(" %10.3e", b[i]);
  printf("\n");
}

//меняет столбцы
static void changestl(double *a, int n, int i, int j)
{
	int q; double c;
	for (q = 0; q < n; q++)
	{
		c = a[q*n+i]; a[q*n+i] = a[q*n+j]; a[q*n+j] = c;
	}
}

//меняет строки начиная с номера start в каждой строке
static void changestr(double *a, int n, int i, int j, int start)
{
	int q; double c; int in = i*n, jn = j*n;
	for (q = start; q < n; q++)
	{
		c = a[in+q]; a[in+q] = a[jn+q]; a[jn+q] = c;
	}
}

//отнимает от i-ой строки j-ую, умноженную на d, начиная с элемента start
static void operation(double *a, int n, int i, int j, double d, int start)
{
	int q; int in = i*n, jn = j*n;
	for (q = start; q < n; q++)
		a[in+q] -= d*a[jn+q];
}


//эта функция с 1 курса, поэтому написана не очень красиво
int get_reversed_matrix_gauss(int n, double *a, double *x, int *r, double eps)
{
	int i,j,imax = 0,jmax = 0,c,kn,in;
	double MAX = 0, aij;
	int k = 0;
	//прямой ход
	while (k < n)
	{
		for (i = k; i < n; i++)
		{
			in = i*n;
			for (j = k; j < n; j++)
			{
				aij = fabs(a[in+j]);
				if (aij > MAX)
				{
					MAX = aij; imax = i; jmax = j;
				}
			}
		}
		if (MAX <= eps)
			return SING_MATRIX;
		if (imax != k)
		{
			changestr(a,n,k,imax,k);
			changestr(x,n,k,imax,0);
		}
		if (jmax != k)
		{
			i = r[k];
			r[k] = r[jmax];
			r [jmax] = i;
			changestl(a,n,k,jmax);
		}
		kn = k*n;
		//проверка что aij не равен нулю была выше
		aij = a[kn+k];
		//делим строки
		for (i = k+1; i < n; i++)
			a[kn+i] /= aij;
		for (i = 0; i < n; i++)
			x[kn+i] /= aij;
		for (i = k+1; i < n; i++)
		{
			aij = a[i*n+k];
			operation(a,n,i,k,aij,k+1);
			operation(x,n,i,k,aij,0);
		}
		MAX = 0;
		k++;
	}
	k = n-1;
	//обратный ход
	while (k > 0)
	{
		for (i = k-1; i >= 0; i--)
		{
			aij = a[i*n+k];
			operation(x,n,i,k,aij,0);
		}
		k--;
	}
	//меняем "столбцы" обратно
	//меняет за квадрат
	
	for (i = 0; i < n; i++)
	{
		c = r[i];
		while (c != i)
		{
			changestr(x,n,c,i,0);
			j = r[c];
			r[c] = r[i];
			r[i] = j;
			c = r[i];
		}
	}
	return SUCCESS;
}

void fill_arranged_items(int *array, int size)
{
	for (int i = 0; i < size; i++)
		array[i] = i;
}


void set_block(double *a, int i, int j, double *b, int n, int m, int height, int width)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			 a[(i*m+h)*n + j*m + w] = b[h*width + w];
		}
	}
	//добавить не для квадратных блоков
}



void get_block(double *a, int i, int j, double *b, int n, int m, int height, int width)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
		  b[h*width + w] = a[(i*m+h)*n + j*m + w];
		}
	}
	//добавить не для квадратных блоков
}

 
void get_vector(double *b, double *res, int ind, int m, int size)
{
	for (int i = 0; i < size; i++)
		res[i] = b[ind*m+i];
}


void set_vector(double *b, double *buf, int ind, int m, int size)
{
	for (int i = 0; i < size; i++)
		b[ind*m + i] = buf[i];
}

void set_unitary_block(double *a, int i, int j, int n, int m, int width, int height)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			a[(i*m+h)*n + j*m + w] = (h == w);
  	}
	}
}

void set_unitary_matrix(double *a, int m)
{
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < m; j++)
	    a[i*m + j] = (i == j);
}



void set_null_block(double *a, int i, int j, int n, int m, int height, int width)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			a[(i*m+h)*n + j*m + w] = 0;
		}
	}
}


int get_reversed_matrix_for_block(double *a, int n, int m, int size, int i, int j, double *res, double eps, double *buf, int * r)
{
	get_block(a, i, j, buf, n, m, size, size);
	fill_arranged_items(r, size);
	set_unitary_matrix(res, size);
	int ret = get_reversed_matrix_gauss(size, buf, res, r, eps); 
	return ret;
}


void change_str_blocks(int *numbers, int i, int j)
{
	int c = numbers[i];
	numbers[i] = numbers[j];
	numbers[j] = c;
}

void arrange_vector_items(double *b, int n, int m, int *numbers)
{
	double *buf = new double[n];
	int s = n/m;
	for (int i = 0; i < s; i++)
	{
		int new_ind = numbers[i];
		for (int j = 0; j < m; j++)
			buf[i*m + j] = b[new_ind*m + j];
	}
	if (n % m != 0)
	{
		int ost = n % m;
		for (int i = 0; i < ost; i++)
			buf[s*m + i] = b[s*m + i];
	}
	copy_vector(b, buf, n);
	delete [] buf;
}

void multiply_block_on_vector(double *buf, double *b, int str_ind, int m, int block_sz, double *b_old)
{
	copy_vector(b_old, b + str_ind*m, block_sz);
	for (int i = 0; i < block_sz; i++)
	{
		b[str_ind*m + i] = 0;
		for (int j = 0; j < block_sz; j++)
		{
			b[str_ind*m + i] += buf[i*block_sz + j]*b_old[j];
		}
	}
}

void multiply_blocks_on_str(double *a, double *mult, int i, int j, int n, int m, double *buf, double *res, double *buf2, double *res2)
{
	int s = n / m;
	int ost = n % m;
	if (i < s)
	{
		for (int w = j; w < s; w++)
		{
			get_block(a, i, w, buf, n, m, m, m);
			//multiply_matrix(mult, buf, res, m, m, m);
			multiply_blocks_m(mult, buf, res, m);
			set_block(a, i, w, res, n, m, m, m);
		}
		if (ost != 0)
		{
			get_block(a, i, s, buf2, n, m, m, ost);
			multiply_matrix(mult, buf2, res2, m, m, ost);
			set_block(a, i, s, res2, n, m, m, ost);
		}
	}
}

void subtract_blocks(double *a, double *b, int n, int m)
{
	for (int i = 0; i < n*m; i++)
		a[i] -= b[i];
}


void subtract_with_multiply(double *a, int from_ind, int ind, int start_stl, double *mult, int n, int m, double *buf, double *res)
{
	int s = n / m;
	int ost = n % m;
	if (from_ind < s)
	{
		for (int w = start_stl; w < s; w++)
		{
			get_block(a, ind, w, buf, n, m, m, m);
			//multiply_matrix(mult, buf, res, m, m, m);
			multiply_blocks_m(mult, buf, res, m);
			get_block(a, from_ind, w, buf, n, m, m, m);
			subtract_blocks(buf, res, m, m);
			set_block(a, from_ind, w, buf, n, m, m, m);
		}
		if (ost != 0)
		{
		  //последний столбец
		  get_block(a, ind, s, buf, n, m, m, ost);
		  multiply_matrix(mult, buf, res, m, m, ost);
		  get_block(a, from_ind, s, buf, n, m, m, ost);
		  subtract_blocks(buf, res, m, ost);
		  set_block(a, from_ind, s, buf, n, m, m, ost); 
		}
	}
	else
	{
		//последняя строка
		for (int w = start_stl; w < s; w++)
		{
		  get_block(a, ind, w, buf, n, m, m, m);
		  multiply_matrix(mult, buf, res, ost, m, m);
		  get_block(a, from_ind, w, buf, n, m, ost, m);
		  subtract_blocks(buf, res, ost, m);
		  set_block(a, from_ind, w, buf, n, m, ost, m);
		}
		//последний маленький блок
		get_block(a, ind, s, buf, n, m, m, ost);
		multiply_matrix(mult, buf, res, ost, m, ost);
		get_block(a, from_ind, s, buf, n, m, ost, ost);
		subtract_blocks(buf, res, ost, ost);
		set_block(a, from_ind, s, buf, n, m, ost, ost);
	}
}


void subtract_with_multiply_vector(double *b, int from_ind, int ind, double *mult, int n, int m, double *buf, double *res, double *res2)
{
	int s = n / m;
	if (from_ind < s)
	{
		get_vector(b, buf, ind, m, m);
		multiply_matrix_and_vector(mult, buf, res, m, m);
		get_vector(b, buf, from_ind, m, m);
		subtract_blocks(buf, res, m, 1);
		set_vector(b, buf, from_ind, m, m);
	}
	else
	{
	  //вычитаем из короткого блока 
	  //умножаем
		int ost = n % m;
		get_vector(b, buf, ind, m, m);
		multiply_matrix_and_vector(mult, buf, res2, ost, m);
		get_vector(b, buf, s, m, ost);
		subtract_blocks(buf, res2, ost, 1);
		set_vector(b, buf, s, m, ost);
	} 
}

int find_block_with_the_min_norm(double *a, int n, int m, int stl, int *numbers, double eps, double *buf, double *buf2, int *r)
{
  int s = n/m;
  double min_norm = 1e308;
  int res = -1;
  for (int i = stl; i < s; i++)
  {
    int ind = numbers[i]; // - индекс реальной матрицы, с которым мы работаем
    int ret = get_reversed_matrix_for_block(a, n, m, m, ind, stl, buf, eps, buf2, r);
    if (ret == SING_MATRIX)
      continue;
    double check = norma_matrix(buf, m, m);
    if (check < min_norm)
    {
      min_norm = check;
      res = i;
    }
  }
  return res;
}


void fill_false(bool *a, int n)
{
  for (int i = 0; i < n; i++)
    a[i] = false;
}

void nill_last_column(double *a, double *b, int n, int m, double *buf, double *res, double *b_last)
{
  int ost = n % m;
  int s = n / m;
  get_vector(b, b_last, s, m, ost);
  for (int i = 0; i < s; i++)
  {
    get_block(a, i, s, buf, n, m, m, ost);
    multiply_matrix_and_vector(buf, b_last, res, m, ost);
    set_null_block(a, i, s, n, m, m, ost);
    get_vector(b, buf, i, m, m);
    subtract_vectors(buf, res, m);
    set_vector(b, buf, i, m, m);
  }
}


static void clean_memory(double *buf, int *numbers, double *buf_mm, double *buf_mm2, double *buf_ostm, double *buf_ostm2, double *buf_ostost, double *buf_ost, double *buf_m, double *buf_m2, int *r_m, int *r_ost)
{
  delete [] buf;
	delete [] numbers;
	delete [] buf_mm;
	delete [] buf_mm2;
	delete [] buf_ostm;
	delete [] buf_ostm2;
	delete [] buf_ostost;
	delete [] buf_ost;
	delete [] buf_m;
	delete [] buf_m2;
	delete [] r_m;
	delete [] r_ost;
}

int block_jordan(double *a, double *b, int n, int m)
{
	int s = n / m;
	int ost = n % m;
	double eps = 1e-16;
	eps *= norma_matrix(a, n, n);
	double *buf = new double[m*m];
	int *numbers = new int[n];
	fill_arranged_items(numbers, n);
	
	//создаём кучу массивов чтобы не создавать их в функциях
	double *buf_mm = new double[m*m];
	double *buf_mm2 = new double[m*m];
	double *buf_ostm = new double[m*ost];
	double *buf_ostm2 = new double[m*ost];
	double *buf_ostost = new double[ost*ost];
	double *buf_ost = new double[ost];
	double *buf_m = new double[m];
	double *buf_m2 = new double[m];
	int *r_m = new int[m];
	int *r_ost = new int[ost];
	//======================================================
	for (int k = 0; k < s; k++)
	{
		int ind_max_numb = find_block_with_the_min_norm(a, n, m, k, numbers, eps, buf_mm, buf_mm2, r_m); //индекс блока в numbers
		if (ind_max_numb < 0) // вырожденная матрица
		{
		  clean_memory(buf, numbers, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
	    return SING_MATRIX;
		}
		int ind_max = numbers[ind_max_numb]; //реальный индекс блока в матрице
		get_reversed_matrix_for_block(a, n, m, m, ind_max, k, buf, eps, buf_mm, r_m);
		set_unitary_block(a, ind_max, k, n, m, m, m);
		multiply_blocks_on_str(a, buf, ind_max, k+1, n, m, buf_mm, buf_mm2, buf_ostm, buf_ostm2);
		multiply_block_on_vector(buf, b, ind_max, m, m, buf_m);
		for (int i = 0; i < s; i++)
		{
			int ind = numbers[i];
			if (ind == ind_max)
				continue;
			get_block(a, ind, k, buf, n, m, m, m);
			if (norma_matrix(buf, m, m) < eps) //если на этой строчке итак элемент в k-ом столбце нулевой
			  continue;
			set_null_block(a, ind, k, n, m, m, m);
			subtract_with_multiply(a, ind, ind_max, k+1, buf, n, m, buf_mm, buf_mm2);
			subtract_with_multiply_vector(b, ind, ind_max, buf, n, m, buf_m, buf_m2, buf_ost);
		}
		if (ind_max_numb != k)
		  change_str_blocks(numbers, ind_max_numb, k);
		
		if (ost != 0)
		{
		  get_block(a, s, k, buf, n, m, ost, m);
			if (norma_matrix(buf, ost, m) < eps) //если на этой строчке итак элемент в k-ом столбце нулевой
			  continue;
			set_null_block(a, s, k, n, m, ost, m);
			subtract_with_multiply(a, s, ind_max, k+1, buf, n, m, buf_mm, buf_mm2);
			subtract_with_multiply_vector(b, s, ind_max, buf, n, m, buf_m, buf_m2, buf_ost);
		}
	}
	if (ost != 0)
	{
	  //поделить эту строчку
	  get_block(a, s, s, buf, n, m, ost, ost);
	  if (norma_matrix(buf, ost, ost) < eps)
	  {
	    clean_memory(buf, numbers, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
	    return SING_MATRIX;
	  }
	  get_reversed_matrix_for_block(a, n, m, ost, s, s, buf, eps, buf_ostost, r_ost);
	  set_unitary_block(a, s, s, n, m, ost, ost);
	  multiply_block_on_vector(buf, b, s, m, ost, buf_ost);
	  nill_last_column(a, b, n, m, buf_ostm, buf_m, buf_ost);
	}
	arrange_vector_items(b, n, m, numbers);
	
	clean_memory(buf, numbers, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
	
	return SUCCESS;
}
