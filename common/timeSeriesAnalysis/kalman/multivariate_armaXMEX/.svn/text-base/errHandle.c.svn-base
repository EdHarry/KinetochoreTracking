#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "carmaFitModel.h"

int errHandle(char *errMessage, int lineNo, char *fileName){
  
  printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

  char errLabel[] = "\n\n*** Error ***\n";

  /* Create string giving location of error */
  char errLoc[100];
  char errLine[] = "Line No.";
  char errFile[] = "in file";
  
  sprintf(errLoc,"%s %d %s %s\n",errLine,lineNo,errFile,fileName);
  

  /* Write error message to std out and to file,
     and append timestamp */

  /* Get the system time and convert to a string */
  time_t timeNow = time(&timeNow); 
  char *timeStampString = ctime(&timeNow);

  /* Print error message and time stamp to std out */
  printf("%sLocation: %s",errLabel,errLoc); 
  printf("%s\n%s\n",errMessage,timeStampString);
  /* Write to file with time */
  int writeErr = 0;

  writeErr += writeLogFile(errLabel);  
  writeErr += writeLogFile(errLoc);
  writeErr += writeLogFile(errMessage);
  writeErr += writeLogFile("\n");
  writeErr += writeLogFile(timeStampString);
  writeErr += writeLogFile("\n");

  if (writeErr != 6){
    printf("Error writing Error to logfile!!\n");
    return 0;
  }

  printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");

  return 1;
}
