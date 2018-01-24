/* 
 * File: prob.c
 * ------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "prob.h"



struct probCDT{
  int arOrd;
  int maOrd;




};



probADT NewProb(int arOrd, int maOrd)
{
  probADT prob;

  prob = (probADT) malloc(sizeof(struct probCDT));
  prob->arOrd = arOrd;
  prob->maOrd = maOrd;

  return prob;
}


void FreeProb(probADT prob)
{
  free(prob);
}
