int *alloc_1d_int_array(int ncolumns)
{
	int *arr = malloc(sizeof(int)* ncolumns);
	return arr;
}

float *alloc_1d_float_array(int ncolumns)
{
	float *arr = malloc(sizeof(float)* ncolumns);
	return arr;
}

double *alloc_1d_double_array(int ncolumns)
{
	double *arr = malloc(sizeof(double)* ncolumns);
	return arr;
}

double **alloc_2darray(int nrows, int ncolumns)
{
	double **arr = malloc(sizeof(double *)* nrows);
	for(int i=0; i < nrows; i++)
		*(arr+i) = malloc(sizeof(double)*ncolumns);
	return arr;
}

int **alloc_2d_int_array(int nrows, int ncolumns)
{
	int **arr = malloc(sizeof(int*)* nrows);
	for(int i=0; i < nrows; i++)
		*(arr+i) = malloc(sizeof(int)*ncolumns);
	return arr;
}

void zero_2darray(double **array, int nrows, int ncolumns)
{
	int i, j;
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < ncolumns; j++)
			array[i][j] = 0;
	}
}

void zero_1darray(double *array, int ncolumns)
{
	for(int i = 0; i < ncolumns; i++)
	 	array[i] = 0;
}


void free2dArray(double **a, int m)
{
	for (int i = 0; i < m; i++)
		free(a[i]);
	free(a);
}

chainList *InitChainList()
{
	chainList *list = malloc(sizeof(*list));
	Element *element = malloc(sizeof(*element));
	if(list == NULL || element == NULL)
		exit(1);

	element->nb = 0;
	element->next = NULL;
	list->first = element;
	list->nbEle =1;

	return list;
}

void InsertEleChainList(chainList *list, int ele)
/* Ajoute un nouvel élément a la première position
 * d'une chainList
 */
{
	Element *nouveau = malloc(sizeof(*nouveau));
	if(nouveau == NULL || list == NULL)
		exit(1);

	nouveau->nb = ele;
	nouveau->next= list->first;
	list->first = nouveau;
	list->nbEle ++;
}

void DelEleChainList(chainList *list)
/* Supprime le premier élément d'une chainList
 */
{
	if (list == NULL)
		exit(1);
	
	if (list->first != NULL)
	{
		Element *del = list->first;
		list->first= list->first->next;
		free(del);
	}
	list->nbEle -= 1 ;
}


void GetChainList(chainList *list)
{
	if(list == NULL)
		exit(1);

	Element *actuel = list->first;
	while(actuel != NULL)
	{
		fprintf(stdout,"%d\n",actuel->nb);
		actuel = actuel->next;
	}
}





