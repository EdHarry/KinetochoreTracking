/*
 * File: prob.c
 * ------------
 * MATLAB input:
 * model.trajOut(1).observations
 * model.trajOut(1).numMissing
 */

/*

struct probCDT{

  double *TRAJ;
  int trajLength;
  int nNodes;
  int numMissing;

  int numParams;
  int arOrderMax;
  int maOrderMax;
  int *topoBIN;
  int *maBIN;

  int nInputs;

  double wnVariance;
};

*/

struct movie {
  double *traj;
  int trajLen;
  int nMissing;
};



/*
 * struct probCDT
 * --------------
 * Data contains an array of movies.
 */

struct probCDT{

  struct movie *data;
  int nMovies;

  int nNodes;
  int nParams;
  int arOrderMax;
  int maOrderMax;
  int *topoBIN;
  int *maBIN;

  double wnVariance;
};




