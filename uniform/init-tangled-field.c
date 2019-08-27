/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author Asger Gronnow (asger.gronnow@usyd.edu.au)
  \date   Nov 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
//#define OUTP_BOV
#define PRS_FROM_TEMP
//#define VELOCITY_MACH
#define TRC_CLOUD TRC
#define TRC_CORE TRC+1
#define TRC_ENV TRC+2
#define NUM_BTERMS 10
#define NORMALIZE_BFIELD
//#define PRS_FROM_FILE
//#define B0_FROM_BETA

#define MICROGAUSS 2.181

/* ********************************************************************* */
double DistanceFromDensity(double HaloDensity, char TableFile[256])
/*!
 * Calculate z distance above galactic midplane for a given halo density
 * based on a table of halo density vs distance, using linear
 * interpolation between data points.
 *
 * \param [in] HaloDensity   Halo baryon number density
 * \param [in] TableFile     Path to file containing table with distance
 *                           in kpc in the first column and halo number
 *                           density per cubic cm in the second column.
 *
 * \return The distance in kpc inferred from the given halo density
 *         or one of the following negative numbers in case of an error:
 *         -1.0 Table file could not be opened
 *         -2.0 Table has too many data points (> 1024)
 *         -3.0 Given density is greater than the maximum value covered
 *          by the table
 *         -4.0 Given density is less than the minimum value covered
 *          by the table
 *********************************************************************** */
{
  static int firstcall = 1;
  static double z_dist [1024];
  static double Density [1024];
  static int Ndata = 0;
  int i = 0;
  if (firstcall)
  {
    FILE* table = fopen(TableFile, "r");
    if (table == NULL) {return -1.0;}
    else
    {
      char line[256];
      double dist;
      double dens;
      double mass;
      double scaled;
      while (fgets(line, sizeof(line), table) != NULL) 
      {
	if (line[0] != '#')
	{
	  sscanf(line, "%lf %lf %lf %lf", &dist, &dens, &mass, &scaled);
	  z_dist[i] = dist;
	  Density[i] = dens;
	  i++;
	}
	if (i > 1023) {return -2.0;}
      }
    }
    fclose(table);
    Ndata = i - 1;
    firstcall = 0;
  }
  if (HaloDensity > Density[0]) {return -3.0;}
  i = 0;
  while (Density[i] > HaloDensity) {i++;}
  if (i > Ndata) {return -4.0;}
  double InterpDistance = z_dist[i] + (z_dist[i+1] - z_dist[i])*(HaloDensity - Density[i])/(Density[i+1] - Density[i]);
  return InterpDistance;
}

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  static int firstcall = 1;
  static int firstcall_prs = 1;

  double CloudRad = g_inputParam[CLOUD_RADIUS];
  double CoreRad = g_inputParam[CORE_RADIUS];
  double CutRad = g_inputParam[CUT_RADIUS];
  double CurRad;
  double TCloud = 1.e3;
  double TWind = 1.e6;
  double muCloud = 1.2;
  double muWind = 1.2;
  double RhoWind = g_inputParam[RHO_W];
  //double RhoCloud = g_inputParam[RHO_C];
  g_gamma = g_inputParam[GAMMA_GAS];
  #ifdef PRS_FROM_TEMP
  double PrsWind = RhoWind*TWind/(KELVIN*muWind);
  #else
  double PrsWind = g_inputParam[PRS_W];
  #endif

  // Cloud pressure given by pressure equilibrium
  double PrsCloud = PrsWind;
  double RhoCloud = 0.1;

  #if PHYSICS == MHD
  unsigned int seed = 3845;
  #ifdef B0_FROM_BETA
  double B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
  #elif BGFIELD != NONE
  // Find distance based on halo density
  // then use extrapolation of linear approximation to the |B|-z relation from Cox (2005) to find background B-field strength
  char densfile[256] = "ISO_MW_halomass.dat";
  double z_dist_kpc = DistanceFromDensity(RhoWind, densfile);
  if (z_dist_kpc < 0) {print1("Error no. %i encountered while attempting to find distance using table %s\n", (int)z_dist_kpc, densfile);}
  double DecayFactor = 0.5*atanh((g_inputParam[BG_B0]-3)/g_inputParam[BG_B0]);
  double B0 = g_inputParam[BG_B0]*(1-tanh(z_dist_kpc*DecayFactor)) * MICROGAUSS;
  if (firstcall) {print1("z_dist: %f, B0: %e, PrsWind: %e, BetaWind: %e\n", z_dist_kpc, B0/MICROGAUSS), 2*PrsWind/(B0*B0);}
  #else
  double B0 = 0;
  #endif
  #ifdef B0_FROM_BETA
  double B0Cloud = sqrt(2*PrsCloud/g_inputParam[BETA_C]);
  #else
  double B0Cloud = g_inputParam[CLOUD_B0] * MICROGAUSS;
  #endif
  double alpha = 10.0/CloudRad;
//  double AvgBfield = 3.07582;//1.14633;
  #endif

  double x = x1 - g_inputParam[CLOUD_CENTREX1];
  double y = x2 - g_inputParam[CLOUD_CENTREX2];
  double z = x3 - g_inputParam[CLOUD_CENTREX3];

  CurRad = EXPAND(x*x , + y*y , + z*z);
  CurRad = sqrt(CurRad);

  v[TRC_CLOUD] = 0;
  v[TRC_CORE] = 0;
  v[TRC_ENV] = 0;

  #if CLOUD_FIELD == FIELD_TANGLED
  int i;

  static double coordvecx[3][NUM_BTERMS];
  static double coordvecy[3][NUM_BTERMS];
  static double coordvecz[3][NUM_BTERMS];

  if (firstcall)
  {
    int j;
    for (i = 0; i < NUM_BTERMS; i++)
    {
      coordvecx[0][i] = 0.5 - rand_r(&seed)/(double)RAND_MAX;
      coordvecy[0][i] = 0.5 - rand_r(&seed)/(double)RAND_MAX;
      coordvecz[0][i] = 0.5 - rand_r(&seed)/(double)RAND_MAX;
      
      coordvecy[1][i] = 0.5 - rand_r(&seed)/(double)RAND_MAX;
      coordvecz[1][i] = 0.5 - rand_r(&seed)/(double)RAND_MAX;
      coordvecx[1][i] = (-coordvecy[0][i]*coordvecy[1][i] - coordvecz[0][i]*coordvecz[1][i])/coordvecx[0][i];

      coordvecx[2][i] = coordvecy[0][i]*coordvecz[1][i] - coordvecz[0][i]*coordvecy[1][i];
      coordvecy[2][i] = coordvecz[0][i]*coordvecx[1][i] - coordvecx[0][i]*coordvecz[1][i];
      coordvecz[2][i] = coordvecx[0][i]*coordvecy[1][i] - coordvecy[0][i]*coordvecx[1][i];

      for (j = 0; j < 3; j++)
      {
	double len = sqrt(coordvecx[j][i]*coordvecx[j][i] + coordvecy[j][i]*coordvecy[j][i] + coordvecz[j][i]*coordvecz[j][i]);
	coordvecx[j][i] /= len;
	coordvecy[j][i] /= len;
	coordvecz[j][i] /= len;
      }
    }
  }
#endif // CLOUD_FIELD == FIELD_TANGLED

  firstcall = 0;
  if (CurRad < CutRad)
  {
    v[RHO] = RhoWind + (RhoCloud - RhoWind)/(1+pow(CurRad/CoreRad,g_inputParam[DENSITY_STEEPNESS]));
    v[PRS] = PrsCloud;
    if (CurRad < CloudRad)
    {
      EXPAND(v[VX1] = g_inputParam[VX1_C]; ,
	     v[VX2] = 0.0; ,
	     v[VX3] = 0.0;)
      v[TRC_CLOUD] = 1.0;
      if (CurRad < CoreRad) {v[TRC_CORE] = 1.0;}
      else {v[TRC_ENV] = 1.0;}
    }
  }
  else
  {
    #ifdef VELOCITY_MACH
    double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsCloud/RhoWind);
    #else
    double VelWind = g_inputParam[VX1_W];
    #endif
    v[RHO] = RhoWind;
    v[PRS] = PrsWind;
    EXPAND(v[VX1] = VelWind; ,
           v[VX2] = 0.0; ,
           v[VX3] = 0.0;)
  }
#if PHYSICS == MHD
  EXPAND(v[AX1] = 0.0; ,
	 v[AX2] = 0.0; ,
	 v[AX3] = 0.0;)
  double Axcloud = 0;
  double Aycloud = 0;
  double Azcloud = 0;

  double CloudBNorm = 1;//0.969310;//0.506901;//0.969902;
#if CLOUD_FIELD == FIELD_TANGLED
  for (i = 0; i < NUM_BTERMS; i++)
  {
    double a = x*coordvecx[0][i] + y*coordvecx[1][i] + z*coordvecx[2][i];
	
    double Bb = sin(alpha*a);
    double Bc = cos(alpha*a);
    double Bx = (coordvecx[1][i]*(Bb*coordvecz[2][i]-Bc*coordvecy[2][i])+coordvecx[2][i]*(Bc*coordvecy[1][i]-Bb*coordvecz[1][i]))/(coordvecx[0][i]*(coordvecy[2][i]*coordvecz[1][i]-coordvecy[1][i]*coordvecz[2][i])+coordvecx[1][i]*(coordvecy[0][i]*coordvecz[2][i]-coordvecy[2][i]*coordvecz[0][i])+coordvecx[2][i]*(coordvecy[1][i]*coordvecz[0][i]-coordvecy[0][i]*coordvecz[1][i]));
    double By = -(coordvecx[0][i]*(Bb*coordvecz[2][i]-Bc*coordvecy[2][i])+coordvecx[2][i]*(Bc*coordvecy[0][i]-Bb*coordvecz[0][i]))/(coordvecx[0][i]*(coordvecy[2][i]*coordvecz[1][i]-coordvecy[1][i]*coordvecz[2][i])+coordvecx[1][i]*(coordvecy[0][i]*coordvecz[2][i]-coordvecy[2][i]*coordvecz[0][i])+coordvecx[2][i]*(coordvecy[1][i]*coordvecz[0][i]-coordvecy[0][i]*coordvecz[1][i]));
    double Bz = (coordvecx[0][i]*(Bb*coordvecz[1][i]-Bc*coordvecy[1][i])+coordvecx[1][i]*(Bc*coordvecy[0][i]-Bb*coordvecz[0][i]))/(coordvecx[0][i]*(coordvecy[2][i]*coordvecz[1][i]-coordvecy[1][i]*coordvecz[2][i])+coordvecx[1][i]*(coordvecy[0][i]*coordvecz[2][i]-coordvecy[2][i]*coordvecz[0][i])+coordvecx[2][i]*(coordvecy[1][i]*coordvecz[0][i]-coordvecy[0][i]*coordvecz[1][i]));

    Axcloud += Bx;
    Aycloud += By;
    Azcloud += Bz;
  }
    
//    Axcloud /= NUM_BTERMS*0.18515/sqrt(2*PrsCloud/g_inputParam[BETA_C]);
//    Aycloud /= NUM_BTERMS*0.18515/sqrt(2*PrsCloud/g_inputParam[BETA_C]);
//    Azcloud /= NUM_BTERMS*0.18515/sqrt(2*PrsCloud/g_inputParam[BETA_C]);

  Axcloud /= alpha;
  Aycloud /= alpha;
  Azcloud /= alpha;
#endif // CLOUD_FIELD == FIELD_TANGLED

  double Axbackground = 0;
  double Aybackground = 0;
  double Azbackground = 0;
#if BG_FIELD == FIELD_TRANSVERSE /* Background field is B=(0, B0, 0) */
  Azbackground = -x*B0;
#if CLOUD_FIELD != FIELD_NONE && CLOUD_FIELD != FIELD_TANGLED
  Azcloud = -x*B0Cloud;
#endif
#elif BG_FIELD == FIELD_PARALLEL /* Background field is B=(B0, 0, 0) */
  Azbackground = y*B0;
#if CLOUD_FIELD != FIELD_NONE && CLOUD_FIELD != FIELD_TANGLED
  Azcloud = y*B0Cloud;
#endif
#elif BG_FIELD == FIELD_OBLIQUE  /* Background field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
  Aybackground = B0*x/sqrt(3);
  Azbackground = B0*(y-x)/sqrt(3);
#endif
/*#if CLOUD_FIELD != FIELD_NONE && CLOUD_FIELD != FIELD_TANGLED
  Aycloud = B0Cloud*x/sqrt(3);
  Azcloud = B0Cloud*(y-x)/sqrt(3);
#endif
*/
  double Wr = 0;
#if CLOUD_FIELD != FIELD_NONE
#if BG_FIELD != FIELD_NONE
  Wr = 1/(1+pow(CurRad/CoreRad,g_inputParam[DENSITY_STEEPNESS]));

  v[AX1] = (Axcloud * Wr + Axbackground * (1 - Wr)) * CloudBNorm;
  v[AX2] = (Aycloud * Wr + Aybackground * (1 - Wr)) * CloudBNorm;
  v[AX3] = (Azcloud * Wr + Azbackground * (1 - Wr)) * CloudBNorm;
#else
  #if CLOUD_FIELD == FIELD_TRANSVERSE /* Cloud field is (0, B0, 0) */
  Axcloud = z*B0Cloud;
  Azcloud = -x*B0Cloud;
  #elif CLOUD_FIELD == FIELD_PARALLEL /* Cloud field is (B0, 0, 0) */
  Azcloud = y*B0Cloud;
  #elif CLOUD_FIELD == FIELD_OBLIQUE  /* Cloud field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
  Aycloud = B0Cloud*x/sqrt(3);
  Azcloud = B0Cloud*(y-x)/sqrt(3);
  #endif

  if (CurRad < CloudRad)
  {
    v[AX1] = Axcloud;
    v[AX2] = Aycloud;
    v[AX3] = Azcloud;
  }
  else
  {
    double r3 = CurRad*CurRad*CurRad;
    double a3 = CloudRad*CloudRad*CloudRad;
    v[AX1] = B0Cloud*(a3/r3)*z;
    v[AX2] = 0;
    v[AX3] = -B0Cloud*(a3/r3)*x;
  }
#endif // BG_FIELD != FIELD_NONE
#endif // CLOUD_FIELD != FIELD_NONE
#endif // PHYSICS == MHD

  #ifdef PRS_FROM_FILE
  if (firstcall_prs)
  {
    int k, input_var[256];
    for (k = 0; k < 256; k++) {input_var[k] = -1;}
    input_var[0] = PRS;
    InputDataSet("grid.out", input_var);
    InputDataRead("prs.init.dbl"," ");
    firstcall_prs = 0;
  }
  InputDataInterpolate(v, x1, x2, x3);
  #endif
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
#ifdef OUTP_BOV
  double dt = 0.1;
  int filenumber = g_time/dt;
  if (prank == 0)
  {
    char* dir;
    char filename[512];
    FILE* BOVout;
    dir = GetOutputDir();
    Output vnames;
    Output* varnames = &vnames;
    varnames->var_name = ARRAY_2D(64,128,char);
    SetDefaultVarNames(varnames);

    int nv = NVAR;
    #ifdef STAGGERED_MHD
     D_EXPAND(
       varnames->var_name[nv]   = "bx1s";
       nv++; ,

       varnames->var_name[nv]   = "bx2s"; 
       nv++; ,

       varnames->var_name[nv]   = "bx3s"; 
       nv++;
     )
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     #if DIMENSIONS == 3   
      varnames->var_name[nv]   = "Ax1";
      nv++;
      varnames->var_name[nv]   = "Ax2";
      nv++;
     #endif
     varnames->var_name[nv]   = "Ax3";
     nv++;
    #endif

    int i;
    for (i = 0; i < nv; i++)
    {
      sprintf(filename, "%s/%s.%04d%s", dir, varnames->var_name[i], filenumber, ".bov");
      BOVout = fopen(filename, "w");
      fprintf(BOVout, "TIME: %f\nDATA_FILE: %s.%04d%s\nDATA_SIZE: %d %d %d\nDATA_FORMAT: DOUBLE\n"
	      "VARIABLE: %s\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n"
	      "BRICK_ORIGIN: %f %f %f\nBRICK_SIZE: %f %f %f",
	      g_time, varnames->var_name[i], filenumber, ".dbl", grid[IDIR].np_int_glob, grid[JDIR].np_int_glob, grid[KDIR].np_int_glob,
	      varnames->var_name[i], g_domBeg[IDIR], g_domBeg[JDIR], g_domBeg[KDIR],
	      g_domEnd[IDIR]-g_domBeg[IDIR], g_domEnd[JDIR]-g_domBeg[JDIR], g_domEnd[KDIR]-g_domBeg[KDIR]);
      fclose(BOVout);
    }

//    filenumber++;
  }
#endif //OUTP_BOV


#ifdef OUTP_CENTER_OF_MASS
  int i, j, k;
    int var = NVAR;
    //  print1("var %f", d->Vc[TRC][middlez][middley][middlex+1]);

    double LocalCenterOfMassX = 0.0;
    double LocalMass = 0.0;
    double GlobalCenterOfMassX = 0.0;
    double GlobalMass = 0.0;
    //double dV = grid[IDIR].dx[0] * grid[JDIR].dx[0] * grid[KDIR].dx[0];

    //    print1("dV %f", grid[IDIR].dx[0] * UNIT_LENGTH);
    //print("rank %i, gbeg %i, gend %i ,ghost %i, globbeg %f globend %f ",prank, grid[IDIR].beg,grid[IDIR].end,grid[IDIR].nghost,grid[IDIR].x_glob[grid[IDIR].beg],grid[IDIR].x_glob[grid[IDIR].end]);
    //    print("rank %i, x1 %f, x2 %f ",prank,(grid[IDIR].x_glob[(grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[0]/2.0), (grid[IDIR].x_glob[IEND + (grid[IDIR].beg - grid[IDIR].nghost)] + grid[IDIR].dx[IEND]/2.0));

    TOT_LOOP(k,j,i)
    {
      LocalMass += d->Vc[RHO][k][j][i];// * dV; //For a uniform grid the volume term is constant cancels out when computing the global center of mass
      int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
      LocalCenterOfMassX += (grid[IDIR].x_glob[GlobalIndex] + grid[IDIR].dx[i]/2.0) * d->Vc[RHO][k][j][i];// * dV;
    }

    /* TODO: MPI STUFF SO THAT THIS WORKS IN PARALLEL
             Hopefully this won't slow the processes down too much. If it does, maybe find a way to only do this ocassionaly rather than at every timestep*/\

    MPI_Allreduce(&LocalCenterOfMassX, &GlobalCenterOfMassX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalMass, &GlobalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    GlobalCenterOfMassX /= GlobalMass;

    double maxrho = -1.0;
    double maxrhoX = 0.0;
    TOT_LOOP(k, j, i)
    {
      if (d->Vc[RHO][k][j][i] > maxrho)
      {
	int GlobalIndex = i + (grid[IDIR].beg - grid[IDIR].nghost);
	maxrho = d->Vc[RHO][k][j][i];
	maxrhoX = grid[IDIR].x_glob[GlobalIndex];
      }
    }

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    double *GlobalmaxrhoX = (double *)malloc(sizeof(double) * numprocs);
    double *Globalmaxrho = (double *)malloc(sizeof(double) * numprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&maxrhoX, 1, MPI_DOUBLE, GlobalmaxrhoX, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&maxrho, 1, MPI_DOUBLE, Globalmaxrho, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    int maxelem = -1.0;
    int elem;
    for (elem=0; elem < numprocs; elem++)
    {
      print("Glob %f ", Globalmaxrho[elem]);
	 if (Globalmaxrho[elem] > Globalmaxrho[maxelem])
	 {
	    maxelem = elem;
	 }
    }

    maxrhoX = GlobalmaxrhoX[maxelem];
    
    free(GlobalmaxrhoX);
    free(Globalmaxrho);

    if (prank == 0)
      {
	char *dir, fname[512];
	FILE *fp;
	dir = GetOutputDir();
	sprintf(fname, "%s/centerofmassx.txt", dir);
	fp = fopen(fname, "a");
	fprintf(fp, "%12.6e  %12.6e %12.6e\n", g_time, GlobalCenterOfMassX, maxrhoX);
	fclose(fp);
      }
#endif //OUTP_CENTER_OF_MASS

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  static int firstcall = 1;
  int i, j, k;

  #if PHYSICS == MHD
  if ((side == 0) && (firstcall))
  {
    firstcall = 0;
    double CutRad = g_inputParam[CLOUD_RADIUS];
    double LocalBmagSum = 0;
    int LocalNCloudCells = 0;

    double LocalBetaSum = 0;
    double GlobalBetaSum = 0;
    int GlobalNCloudCells = 0;

    DOM_LOOP(k,j,i)
    {
      int GlobalIndex_i = i + (grid[IDIR].beg - grid[IDIR].nghost);
      int GlobalIndex_j = j + (grid[JDIR].beg - grid[JDIR].nghost);
      int GlobalIndex_k = k + (grid[KDIR].beg - grid[KDIR].nghost);
      
      double x = grid[IDIR].x_glob[GlobalIndex_i] - g_inputParam[CLOUD_CENTREX1];
      double y = grid[JDIR].x_glob[GlobalIndex_j] - g_inputParam[CLOUD_CENTREX2];
      double z = grid[KDIR].x_glob[GlobalIndex_k] - g_inputParam[CLOUD_CENTREX3];

      double CurRad = EXPAND(x*x , + y*y , + z*z);
      CurRad = sqrt(CurRad);

      if (CurRad < CutRad)
      {
	LocalNCloudCells++;
	//LocalBmagSum += sqrt(d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i] + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i] + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]);
	#ifdef B0_FROM_BETA
	LocalBetaSum += 2*d->Vc[PRS][k][j][i]/(d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i] + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i] + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]);
	#else
	LocalBetaSum += sqrt(d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i] + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i] + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]);
	#endif
//	printf("Term: %f\n",2*d->Vc[PRS][k][j][i]/(d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i] + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i] + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]));
      }
    }

    //MPI_Allreduce(&LocalBmagSum, &GlobalBmagSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalBetaSum, &GlobalBetaSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&LocalNCloudCells, &GlobalNCloudCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //double AvgBmag = GlobalBmagSum/GlobalNCloudCells;
    #ifdef B0_FROM_BETA
    double NormConst  = sqrt((GlobalBetaSum/GlobalNCloudCells)/g_inputParam[BETA_C]);
    #else
    double NormConst = g_inputParam[CLOUD_B0]*MICROGAUSS/(GlobalBetaSum/GlobalNCloudCells);
    #endif
    //printf("AvgBmag %f, NCells %i, BmagSum %f, Norm. const. %f\n", AvgBmag, GlobalNCloudCells, GlobalBmagSum, sqrt(2*g_inputParam[PRS_C]/g_inputParam[BETA_C])/AvgBmag);
    print1("AvgBeta %f, NCells %i, BetaSum %f, Norm. const. %f\n", GlobalBetaSum/GlobalNCloudCells, GlobalNCloudCells, GlobalBetaSum, NormConst);

    #ifdef NORMALIZE_BFIELD
    CutRad = g_inputParam[CUT_RADIUS];
    DOM_LOOP(k,j,i)
    {
      int GlobalIndex_i = i + (grid[IDIR].beg - grid[IDIR].nghost);
      int GlobalIndex_j = j + (grid[JDIR].beg - grid[JDIR].nghost);
      int GlobalIndex_k = k + (grid[KDIR].beg - grid[KDIR].nghost);
      
      double x = grid[IDIR].x_glob[GlobalIndex_i] - g_inputParam[CLOUD_CENTREX1];
      double y = grid[JDIR].x_glob[GlobalIndex_j] - g_inputParam[CLOUD_CENTREX2];
      double z = grid[KDIR].x_glob[GlobalIndex_k] - g_inputParam[CLOUD_CENTREX3];

      double CurRad = EXPAND(x*x , + y*y , + z*z);
      CurRad = sqrt(CurRad);

      #ifdef PRS_FROM_TEMP
      double TWind = 1.e6;
      double muWind = 1.2;
      double RhoWind = g_inputParam[RHO_W];
      double PrsCloud = RhoWind*TWind/(KELVIN*muWind); //Same pressure as for wind because of pressure equilibrium
      #else
      double PrsCloud = g_inputParam[PRS_C];
      #endif

      if (CurRad < CutRad)
      {
	d->Vc[BX1][k][j][i] *= NormConst;//sqrt(2*PrsCloud/g_inputParam[BETA_C])/AvgBmag;
	d->Vc[BX2][k][j][i] *= NormConst;//sqrt(2*PrsCloud/g_inputParam[BETA_C])/AvgBmag;
	d->Vc[BX3][k][j][i] *= NormConst;//sqrt(2*PrsCloud/g_inputParam[BETA_C])/AvgBmag;
	d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY; /* Tell PLUTO to update cell. MUST be present. */
      }
    }
    #endif //NORMALIZE_BFIELD
  }
#endif //PHYSICS == MHD

  // set wind (inflow) properties
  if (side == X1_BEG)
  {
    if (box->vpos == CENTER)
    {
      double TWind = 1.e6;
      double muWind = 1.2;
      double RhoWind = g_inputParam[RHO_W];
      double RhoCloud = g_inputParam[RHO_C];
      #ifdef PRS_FROM_TEMP
      double PrsWind = RhoWind*TWind/(KELVIN*muWind);
      #else
      double PrsWind = g_inputParam[PRS_W];
      #endif
      #ifdef VELOCITY_MACH
      double VelWind = g_inputParam[VX1_W] * sqrt(g_gamma*PrsWind/RhoWind);
      #else
      double VelWind = g_inputParam[VX1_W];
      #endif

      BOX_LOOP(box,k,j,i)
      {
	d->Vc[RHO][k][j][i] = RhoWind;
	double PrsBCloud = 0;
	d->Vc[PRS][k][j][i] = PrsWind;
	EXPAND(d->Vc[VX1][k][j][i] = VelWind; ,
	  d->Vc[VX2][k][j][i] = 0.0; ,
	  d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[TRC_CLOUD][k][j][i] = 0.0;
        d->Vc[TRC_CORE][k][j][i] = 0.0;
        d->Vc[TRC_ENV][k][j][i] = 0.0;
#if PHYSICS == MHD
	#ifdef B0_FROM_BETA
	dobule B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
	#else
	double B0 = g_inputParam[BG_B0]*MICROGAUSS;
	#endif
	d->Vc[BX1][k][j][i] = 0;
	d->Vc[BX2][k][j][i] = 0;
	d->Vc[BX3][k][j][i] = 0;
#if BG_FIELD == FIELD_TRANSVERSE /* Background field is B=(0, B0, 0) */
	d->Vc[BX2][k][j][i] = B0;
#elif BG_FIELD == FIELD_PARALLEL /* Background field is B=(B0, 0, 0) */
	d->Vc[BX1][k][j][i] = B0;
#elif BG_FIELD == FIELD_OBLIQUE  /* Background field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
	d->Vc[BX1][k][j][i] = B0/sqrt(3);
	d->Vc[BX2][k][j][i] = B0/sqrt(3);
	d->Vc[BX3][k][j][i] = B0/sqrt(3);
#endif
#endif // PHYSICS == MHD
      }
    }
#ifdef STAGGERED_MHD
    else if ((box->vpos == X1FACE) || (box->vpos == X2FACE) || (box->vpos == X3FACE))
    {
      double TWind = 1.e6;
      double muWind = 1.2;
      double RhoWind = g_inputParam[RHO_W];
      #ifdef PRS_FROM_TEMP
      double PrsWind = RhoWind*TWind/(KELVIN*muWind);
      #else
      double PrsWind = g_inputParam[PRS_W];
      #endif
      #ifdef B0_FROM_BETA
      double B0 = sqrt(2*PrsWind/g_inputParam[BETA_W]);
      #else
      double B0 = g_inputParam[BG_B0]*MICROGAUSS;
      #endif
      BOX_LOOP(box,k,j,i)
      {
	if (box->vpos == X1FACE) {d->Vs[BX1s][k][j][i] = 0;}
	else if (box->vpos == X2FACE) {d->Vs[BX2s][k][j][i] = 0;}
	else if (box->vpos == X3FACE) {d->Vs[BX3s][k][j][i] = 0;}
#if BG_FIELD == FIELD_TRANSVERSE /* Background field is B=(0, B0, 0) */
	if (box->vpos == X2FACE) {d->Vs[BX2s][k][j][i] = B0;}
#elif BG_FIELD == FIELD_PARALLEL /* Background field is B=(B0, 0, 0) */
	if (box->vpos == X1FACE) {d->Vs[BX1s][k][j][i] = B0;}
#elif BG_FIELD == FIELD_OBLIQUE  /* Background field is B=(B0/sqrt(3), B0/sqrt(3), B0/sqrt(3)) */
	if (box->vpos == X1FACE) {d->Vs[BX1s][k][j][i] = B0/sqrt(3);}
        if (box->vpos == X2FACE) {d->Vs[BX2s][k][j][i] = B0/sqrt(3);}
        if (box->vpos == X3FACE) {d->Vs[BX3s][k][j][i] = B0/sqrt(3);}
#endif
      }
    }
    
#endif // STAGGERED_MHD
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
