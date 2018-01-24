#include <stdlib.h>

void createSimplex(double *paramV, int numParams, double edgeLength, double **p)
{
  int i, j;

  for (i = 1; i <= numParams + 1; i++){
    for (j = 1; j <= numParams; j++){

      if (i == j + 1)
	p[i][j] = paramV[j] + edgeLength;
      else
	p[i][j] = paramV[j];

    }
  }
}
