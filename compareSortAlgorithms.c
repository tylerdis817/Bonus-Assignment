#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int extraMemoryAllocated;


//Below 3 functions were provided
void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}



void heapHelp(int arr[], int n, int i) {
    int largest = i; //Initialize largest as root
    int l = 2*i + 1; //left
    int r = 2*i + 2; //right

    //If left child is larger than root
    if (l < n && arr[l] > arr[largest])
        largest = l;

    //If right child is larger than largest so far
    if (r < n && arr[r] > arr[largest])
        largest = r;

    //Checks largest is not the root
    if (largest != i) {
        // Swap
        int swap = arr[i];
        arr[i] = arr[largest];
        arr[largest] = swap;

        //Recursion heap the affected sub-tree
        heapHelp(arr, n, largest);
    }
}


// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapSort(int arr[], int n)
{
     //Build heap
    for (int i = n / 2 - 1; i >= 0; i--)
        heapHelp(arr, n, i);

    //One by one extract an element
    for (int i=n-1; i>0; i--) {
        //Move current root to end
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        //Call max heapHelp on the reduced heap
        heapHelp(arr, i, 0);
    }
}


// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
	if (l < r) {
        //Find the middle
        int m = l + (r - l) / 2;

        //Sort first and second halves
        mergeSort(pData, l, m);
        mergeSort(pData, m + 1, r);

        //Merge the sorted halves
        int n1 = m - l + 1;
        int n2 = r - m;

        //Set memory for the temporary arrs
        int *L = (int *)Alloc(n1 * sizeof(int));
        int *R = (int *)Alloc(n2 * sizeof(int));

        //Cpoies data to temp arrays L and R
        for (int i = 0; i < n1; i++)
            L[i] = pData[l + i];
        for (int j = 0; j < n2; j++)
            R[j] = pData[m + 1 + j];

        //Merges the temporary arrays back into arr
        int i = 0, j = 0, k = l;
        while (i < n1 && j < n2) {
            if (L[i] <= R[j]) {
                pData[k] = L[i];
                i++;
            } else {
                pData[k] = R[j];
                j++;
            }
            k++;
        }

        //Copy any remaining elements of L then R
        while (i < n1)
		{

            pData[k] = L[i];
            i++;
            k++;
			extraMemoryAllocated += sizeof(int);
        }
        while (j < n2)
		{
            pData[k] = R[j];
            j++;
            k++;
			extraMemoryAllocated += sizeof(int);

        }

        //Free the allocated memory
        DeAlloc(L);
        DeAlloc(R);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{

    int i, key, j;
	//Starts from second element in array
    for (i = 1; i < n; i++)
    {

        key = pData[i];
        j = i - 1;

		//Moves larger elements than key a position ahead
        while (j >= 0 && pData[j] > key)
        {
            pData[j + 1] = pData[j];
            j = j - 1;
			extraMemoryAllocated += sizeof(int);
        }
		//Puts key in correct spot
        pData[j + 1] = key;
    }


}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{

   int i, j;
   for (i = 0; i < n-1; i++)       
       for (j = 0; j < n-i-1; j++) 
           if (pData[j] > pData[j+1])
           {

              int temp = pData[j];
              pData[j] = pData[j+1];
              pData[j+1] = temp;
			  extraMemoryAllocated += sizeof(int);
           }
}



// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
for (int i = 0; i < n - 1; i++) {
        // Find the minimum element in unsorted array
        int minIdx = i;
        for (int j = i + 1; j < n; j++)
            if (pData[j] < pData[minIdx])
                minIdx = j;

        //Swaps the smallest element with the first element of the unsorted portion
        if (minIdx != i) {
            int temp = pData[minIdx];
            pData[minIdx] = pData[i];
            pData[i] = temp;
			extraMemoryAllocated += sizeof(int);
        }
    }
}



// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	//Opens files
	FILE* inFile = fopen(inputFileName, "r");
    int dataSz = 0;
	//Gets data from file
    if (inFile) {
        fscanf(inFile, "%d", &dataSz);
        *ppData = (int *)Alloc(sizeof(int) * dataSz);
        for (int i = 0; i < dataSz; i++) {
            fscanf(inFile, "%d", (*ppData) + i);
        }
		//Closes file
        fclose(inFile);
    }
    return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

//main was provided
int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz-1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
