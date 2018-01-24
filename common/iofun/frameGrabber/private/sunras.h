/* general type definitions for the mex{Read/Write}SunRaster functions */

#ifndef __SUNRAS__
#define __SUNRAS__

int readImage(unsigned char **,const char *,int *,int *);

#define RAS_MAGIC       0x59a66a95
#define RT_STANDARD     1       /* Raw pixrect image in 68000 byte order */
                                /* Sun supported ras_maptype's */
#define RMT_NONE        0       /* ras_maplength is expected to be 0 */
#define RMT_EQUAL_RGB   1       /* red[ras_maplength/3],green[],blue[] */


/* error codes */
#define FILEOPENFAILED -1
#define UNSUPPORTEDDEPTH -2
#define NOTSTANDARDSUNRAS -4

#endif
