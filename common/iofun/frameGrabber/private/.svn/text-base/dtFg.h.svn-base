/* general type definitions for the mex{Read/Write}SunRaster functions */

#ifndef __DTFG__
#define __DTFG__

/* header files required for the DataTranslation framegrabber */
#include <windows.h>     /* required because DT lib needs Windows macros */
#include <windowsx.h>

#include "olwintyp.h"
#include "olimgapi.h"
#include "olfgapi.h"


/* types */
typedef struct{
  char alias[OLC_MAX_ALIAS_STR_SIZE];
  OLT_IMG_DEV_ID devId;
  char devName[OLC_MAX_DEVICE_NAME_STR_SIZE + 1];
  ULNG nInpSrc;
  USHRT inpSrc;
  ULNG  fmTypeCaps;
  ULNG  pixDepth;
  ULNG  fmHeight, fmWidth , fmType;
  OLT_FG_FRAME_ID fmId;
  int isOpen;
}DEVINFO;


/* function declarations */
int openDTFg(int*);
int closeDTFg(int*);
int selectDTFg(LPOLT_IMGDEVINFO,int);

int getCtrlDTFg(int*);
void printCtrlDTFg(void);
void printFrameInfoDTFg(LPOLT_FG_FRAME_INFO);
void printImgDevInfoDTFg(LPOLT_IMGDEVINFO);

int grabFromDTFg(unsigned char*,int[4],ULNG,int*);


/* constants */
#define MSGBUFSIZE 128

/* error codes */
#define DTSPECIAL -1
#define GETDEVICECOUNTFAILED -2
#define NODTDEVICEFOUND -3
#define ROIDEFERROR -4

/* errors in Win32 global memory managment 
   (see GlobalAlloc(), GlobalFree(), etc.) */
#define ALLOCATIONERROR -1001
#define LOCKERROR -1002

#endif
