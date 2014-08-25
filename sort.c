#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Simple shell sort, from http://rosettacode.org/wiki/Sorting_algorithms/Shell_sort
void sort_func(float *a, int *jlist, int n)
{
  int h, i, j, jtmp;
  float atmp;
  for (h = n; h /= 2;) {
    for (i = h; i < n; i++) {
      atmp = a[i];
      jtmp = jlist[i];
      for (j = i; j >= h && atmp < a[j - h]; j -= h) {
	a[j] = a[j - h];
	jlist[j] = jlist[j - h];
      }
      a[j] = atmp;
      jlist[j] = jtmp;
    }
  }
}
 
