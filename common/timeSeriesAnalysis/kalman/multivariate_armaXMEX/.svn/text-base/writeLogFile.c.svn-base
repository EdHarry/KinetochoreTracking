#include <stdlib.h>
#include <stdio.h>
#include "carmaFitModel.h"

int writeLogFile(char *logMessage){

  FILE *logFile;
  logFile = fopen("noLABtestLogFile.txt" ,"a");
 
  if (  !fputs(logMessage,logFile) ){
    printf("writeLogFile: Cannot open log file!\n");
    return 0;
  }
  
  if ( fclose(logFile) == EOF ){
    printf("writeLogFile: Cannot close log file!\n");
    return 0;
  }
  return 1;
}
