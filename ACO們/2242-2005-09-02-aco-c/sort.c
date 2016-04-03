

#include<stdio.h>

void swap(long int v[], long int i, long int j);
void sort(long int v[], long int left, long int right);


void swap(long int v[], long int i, long int j)
/*    
      FUNCTION:       auxiliary routine for sorting an integer array  
      INPUT:          array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of array are swapped
*/
{
  long int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}




void sort(long int v[], long int left, long int right)
/*    
      FUNCTION:       recursive routine (quicksort) for sorting an array  
      INPUT:          one array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  long int k, last;

  if (left >= right) 
    return;
  swap(v, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap(v, ++last, k);
  swap(v, left, last);
  sort(v, left, last);
  sort(v, last+1, right);
}


