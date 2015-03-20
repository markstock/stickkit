/*
 * stickkit.c
 * bob.c
 *
 * copyright 2014 Mark J. Stock mstock@umich.edu
 *
 * a program to read, manipulate, and write string, segmented, or
 * graph-like data in fewer than 16 space dimensions
 *
 * compile with
 *   make
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "stickkit.h"

#define MIN(a,b) ({ \
    typeof(a) _a_temp_; \
    typeof(b) _b_temp_; \
    _a_temp_ = (a); \
    _b_temp_ = (b); \
    _a_temp_ = _a_temp_ > _b_temp_ ? _b_temp_ : _a_temp_; \
    })

#define MAX(a,b) ({ \
    typeof(a) _a_temp_; \
    typeof(b) _b_temp_; \
    _a_temp_ = (a); \
    _b_temp_ = (b); \
    _a_temp_ = _a_temp_ < _b_temp_ ? _b_temp_ : _a_temp_; \
    })

/*Function to find minimum of x and y*/
int min(int x, int y)
{
  return y ^ ((x ^ y) & -(x < y));
}
 
/*Function to find maximum of x and y*/
int max(int x, int y)
{
  return x ^ ((x ^ y) & -(x < y));
}

// needed externs
extern int update_node_group_stats (node_group_ptr, int);

// file-local declarations
unsigned char*** allocate_3d_array_b(int nx, int ny, int nz) {

   // allocate an array of nx pointers, one for each plane
   unsigned char ***array = (unsigned char ***)malloc(nx * sizeof(char **));

   // are we dealing with more than 2B bytes?
   // On x86-64 Linux, INT_MAX is 2147483647, LONG_MAX is 9223372036854775807
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(char);

   if (totBytes > (long int)INT_MAX/2) {
      // all data is larger than 2GB, allocate it in planes
      fprintf(stderr,"  voxel array is plane-by-plane\n");

      for (int i=0; i<nx; i++) {
         array[i] = (unsigned char **)malloc(ny * sizeof(char *));
         array[i][0] = (unsigned char *)malloc(ny * nz * sizeof(char));
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }

   } else {
      // we can fit it all in one malloc (it's under 2GB)
      fprintf(stderr,"  voxel array is monolithic\n");

      array[0] = (unsigned char **)malloc(nx * ny * sizeof(char *));
      for (int i=1; i<nx; i++)
         array[i] = array[0] + i * ny;

      array[0][0] = (unsigned char *)malloc(nx * ny * nz * sizeof(char));
      for (int i=0; i<nx; i++) {
         if (i!=0)
            array[i][0] = array[0][0] + i * ny * nz;
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }
   }

   return(array);
}

int free_3d_array_b(unsigned char*** array, int nx, int ny, int nz){
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(char);
   if (totBytes > (long int)INT_MAX/2) {
      for (int i=0; i<nx; i++) {
         free(array[i][0]);
         free(array[i]);
      }
   } else {
      free(array[0][0]);
      free(array[0]);
   }
   free(array);
   return(0);
}

unsigned int*** allocate_3d_array_i(int nx, int ny, int nz) {

   // allocate an array of nx pointers, one for each plane
   unsigned int ***array = (unsigned int ***)malloc(nx * sizeof(int **));

   // are we dealing with more than 2B bytes?
   // On x86-64 Linux, INT_MAX is 2147483647, LONG_MAX is 9223372036854775807
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(int);

   if (totBytes > (long int)INT_MAX/2) {
      // all data is larger than 2GB, allocate it in planes
      fprintf(stderr,"  voxel array is plane-by-plane\n");

      for (int i=0; i<nx; i++) {
         array[i] = (unsigned int **)malloc(ny * sizeof(int *));
         array[i][0] = (unsigned int *)malloc(ny * nz * sizeof(int));
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }

   } else {
      // we can fit it all in one malloc (it's under 2GB)
      fprintf(stderr,"  voxel array is monolithic\n");

      array[0] = (unsigned int **)malloc(nx * ny * sizeof(int *));
      for (int i=1; i<nx; i++)
         array[i] = array[0] + i * ny;

      array[0][0] = (unsigned int *)malloc(nx * ny * nz * sizeof(int));
      for (int i=0; i<nx; i++) {
         if (i!=0)
            array[i][0] = array[0][0] + i * ny * nz;
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }
   }

   return(array);
}

int free_3d_array_i(unsigned int*** array, int nx, int ny, int nz){
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(int);
   if (totBytes > (long int)INT_MAX/2) {
      for (int i=0; i<nx; i++) {
         free(array[i][0]);
         free(array[i]);
      }
   } else {
      free(array[0][0]);
      free(array[0]);
   }
   free(array);
   return(0);
}

/*
 * Write a 3D brick of bytes
 */
int write_bob_file_from_uchar(FILE* ofp, unsigned char*** z, int nx, int ny, int nz) {

   /* write header */
   fwrite(&nx,sizeof(int),1,ofp);
   fwrite(&ny,sizeof(int),1,ofp);
   fwrite(&nz,sizeof(int),1,ofp);

   /* write the data */
   for (int i=0; i<nx; i++) {
      fwrite(&z[i][0][0],sizeof(unsigned char),ny*nz,ofp);
   }

   /* return 0 if all went well */
   return(0);
}

//
// min dist code from
// http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
//
double minimum_distance(double vx, double vy, double vz,
                        double wx, double wy, double wz,
                        double px, double py, double pz) {
  //fprintf(stderr,"v %g %g %g\n",vx,vy,vz);
  //fprintf(stderr,"w %g %g %g\n",wx,wy,wz);
  //fprintf(stderr,"p %g %g %g\n",px,py,pz);
  // Return minimum distance between line segment vw and point p
  // i.e. |w-v|^2 -  avoid a sqrt
  const double l2 = pow(vx-wx,2) + pow(vy-wy,2) + pow(vz-wz,2);
  // v == w case
  if (l2 == 0.0) return sqrt( pow(vx-px,2) + pow(vy-py,2) + pow(vz-pz,2) );
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  const double t = ( (px-vx)*(wx-vx) + (py-vy)*(wy-vy) + (pz-vz)*(wz-vz) ) / l2;
  // Beyond the 'v' end of the segment
  if (t < 0.0) return sqrt( pow(vx-px,2) + pow(vy-py,2) + pow(vz-pz,2) );
  // Beyond the 'w' end of the segment
  else if (t > 1.0) return sqrt( pow(wx-px,2) + pow(wy-py,2) + pow(wz-pz,2) );
  // Projection falls on the segment
  const double jx = vx + t * (wx - vx);
  const double jy = vy + t * (wy - vy);
  const double jz = vz + t * (wz - vz);
  return sqrt( pow(jx-px,2) + pow(jy-py,2) + pow(jz-pz,2) );
}

double min_dist_with_rad(double vx, double vy, double vz, double vr,
                         double wx, double wy, double wz, double wr,
                         double px, double py, double pz) {
  // Return minimum distance between fat line segment vw and point p
  // i.e. |w-v|^2 -  avoid a sqrt
  const double l2 = pow(vx-wx,2) + pow(vy-wy,2) + pow(vz-wz,2);
  const double dl = sqrt(l2);
  const double dvp = sqrt( pow(vx-px,2) + pow(vy-py,2) + pow(vz-pz,2) );
  const double dwp = sqrt( pow(wx-px,2) + pow(wy-py,2) + pow(wz-pz,2) );
  // old: v == w case
  //if (l2 == 0.0) return dvp - fmax(vr, wr);
  // new: one radius overwhelms the entire segment and other endcap
  if (vr > dl+wr) return dvp - vr;
  if (wr > dl+vr) return dwp - wr;
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  double t = ( (px-vx)*(wx-vx) + (py-vy)*(wy-vy) + (pz-vz)*(wz-vz) ) / l2;
  //fprintf(stdout,"\nold t %g\n",t);
  if (fabs(vr-wr) / dl > 1.e-5) {
    // segment is enough of a cone to matter
    // must shift t to account for cone
    // project to centerline, regardless of t
    const double jx = vx + t * (wx - vx);
    const double jy = vy + t * (wy - vy);
    const double jz = vz + t * (wz - vz);
    const double pt = sqrt( pow(jx-px,2) + pow(jy-py,2) + pow(jz-pz,2) );
    // push t one direction or another
    t += (pt/dl) * (wr-vr) / sqrt(l2 - pow(wr-vr,2));
    //fprintf(stdout,"new t %g\n",t);
  }
  //fprintf(stdout,"dvp %g\n",dvp);
  //fprintf(stdout,"dwp %g\n",dwp);
  // Beyond the 'v' end of the segment
  if (t < 0.0) return dvp - vr;
  // Beyond the 'w' end of the segment
  else if (t > 1.0) return dwp - wr;
  // Projection falls on the segment
  const double jx = vx + t * (wx - vx);
  const double jy = vy + t * (wy - vy);
  const double jz = vz + t * (wz - vz);
  // and find the local radius
  const double thisRad = vr + t*(wr-vr);
  return sqrt( pow(jx-px,2) + pow(jy-py,2) + pow(jz-pz,2) ) - thisRad;
}

//
// rasterize all segments to a regular 3D grid of bytes
//
int write_bob(FILE* ofp, seg_group_ptr thisSG, double dx) {

  int nx, ny, nz;
  double start[3];
  double size[3];
  unsigned char*** dat = NULL;

  // wrong if only 2D
  if (thisSG->dim != 3) {
    fprintf(stderr,"Will not write 2D segments to brick-of-bytes.\n");
    fflush(stderr);
    return(1);
  }

  // how big is this domain?
  //(void) update_node_group_stats(thisSG->nodes, thisSG->dim);
  fprintf(stderr,"  node minima");
  for (int i=0; i<3; i++) fprintf(stderr," %g",thisSG->nodes->min[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  node maxima");
  for (int i=0; i<3; i++) fprintf(stderr," %g",thisSG->nodes->max[i]);
  fprintf(stderr,"\n");

  // set build domain 10% larger (will this be enough to capture wide lines?)
  for (int i=0; i<3; i++) size[i] = thisSG->nodes->max[i] - thisSG->nodes->min[i];
  const double maxSize = fmax(fmax(size[0],size[1]),size[2]);
  double maxRad = thisSG->radii->max;
  if (maxRad < 0.0) maxRad = thisSG->radius;
  if (maxRad < 0.0) maxRad = maxSize / 100.0;
  maxRad += 2.0*dx;
  for (int i=0; i<3; i++) start[i] = thisSG->nodes->min[i] - maxRad;
  for (int i=0; i<3; i++) size[i] += 2.0*maxRad;

  // compute reasonable bounds
  fprintf(stderr,"  cell size %g\n",dx);
  nx = size[0] / dx;
  ny = size[1] / dx;
  nz = size[2] / dx;
  fprintf(stderr,"  brick will be %d x %d x %d\n",nx,ny,nz);

  // sanity check on bob size
  if (nx*(float)ny*nz > 1.e+11 || nx > 100000 || ny > 100000 || nz > 100000) {
    fprintf(stderr,"Are you sure you want to write a brick-of-bytes file that large?\n");
    fprintf(stderr,"  Hit enter to continue, control-c to quit.\n");
    fflush(stderr);
    char buf[2];
    fread(&buf, sizeof(char), 1, stdin);
  }

  // allocate space for the brick
  dat = allocate_3d_array_b(nx, ny, nz);

  // zero out
  for (int i=0; i<nx; i++) {
  for (int j=0; j<ny; j++) {
  for (int k=0; k<nz; k++) {
    dat[i][j][k] = 0;
  }
  }
  }

  fprintf(stderr,"  computing");
  fflush(stderr);

  if (FALSE) {
    // run some tests
    double val = min_dist_with_rad(1.0,0.0,0.0,1.0, 3.0,0.0,0.0,2.0, 1.0,2.0,0.0);
    fprintf(stdout,"dist to 1,2 is %g\n",val);
    val = min_dist_with_rad(1.0,0.0,0.0,1.0, 3.0,0.0,0.0,2.0, 2.0,2.0,0.0);
    fprintf(stdout,"dist to 2,2 is %g\n",val);
    val = min_dist_with_rad(1.0,0.0,0.0,1.0, 3.0,0.0,0.0,2.0, 3.0,2.0,0.0);
    fprintf(stdout,"dist to 3,2 is %g\n",val);

    val = min_dist_with_rad(1.0,0.0,0.0,2.0, 3.0,0.0,0.0,1.0, 1.0,2.0,0.0);
    fprintf(stdout,"dist to 1,2 is %g\n",val);
    val = min_dist_with_rad(1.0,0.0,0.0,2.0, 3.0,0.0,0.0,1.0, 2.0,2.0,0.0);
    fprintf(stdout,"dist to 2,2 is %g\n",val);
    val = min_dist_with_rad(1.0,0.0,0.0,2.0, 3.0,0.0,0.0,1.0, 3.0,2.0,0.0);
    fprintf(stdout,"dist to 3,2 is %g\n",val);
    exit(0);
  }

  // then, rasterize all segments
  int nlines = 0;
  seg_ptr curr = thisSG->first;
  while (curr) {

    double rad1;
    if (curr->r[0]) {
      rad1 = curr->r[0]->r;
    } else {
      rad1 = thisSG->radius;
    }
    // insure that we have a positive and appropriate radius
    if (rad1 < 0.0) rad1 = 2.0*dx;
    rad1 /= dx;

    double rad2;
    if (curr->r[1]) {
      rad2 = curr->r[1]->r;
    } else {
      rad2 = thisSG->radius;
    }
    // insure that we have a positive and appropriate radius
    if (rad2 < 0.0) rad2 = 2.0*dx;
    rad2 /= dx;

    // scale the segment into grid coords
    const double x1 = (curr->n[0]->x[0] - start[0]) / dx;
    const double y1 = (curr->n[0]->x[1] - start[1]) / dx;
    const double z1 = (curr->n[0]->x[2] - start[2]) / dx;
    const double x2 = (curr->n[1]->x[0] - start[0]) / dx;
    const double y2 = (curr->n[1]->x[1] - start[1]) / dx;
    const double z2 = (curr->n[1]->x[2] - start[2]) / dx;

    // find x,y,z range affected by this segment
    const int imin = max((int)floor(fmin(x1-rad1, x2-rad2)) - 1, 0);
    const int imax = min((int)ceil(fmax(x1+rad1, x2+rad2)) + 1, nx);
    const int jmin = max((int)floor(fmin(y1-rad1, y2-rad2)) - 1, 0);
    const int jmax = min((int)ceil(fmax(y1+rad1, y2+rad2)) + 1, ny);
    const int kmin = max((int)floor(fmin(z1-rad1, z2-rad2)) - 1, 0);
    const int kmax = min((int)ceil(fmax(z1+rad1, z2+rad2)) + 1, nz);
    //printf("kmax  %d %d %d\n",(int)ceil(fmax(z1+rad1, z2+rad2)), nz, kmax);

    // loop over that subblock
    for (int i=imin; i<imax; i++) {
    for (int j=jmin; j<jmax; j++) {
    for (int k=kmin; k<kmax; k++) {
      // how far is this node from the segment, in voxels?
      double thisDist = 2.0;

      if (FALSE) {
        // do it the old way, assuming constant radius
        thisDist = minimum_distance(x1,y1,z1, x2,y2,z2, (double)i+0.5,(double)j+0.5,(double)k+0.5) - rad1;
        // set the scaled char value: 255=inside, 0=more than 1 voxel away from surface
        //if (thisDist - rad1 < -1.0) thisChar = 255;
        //else if (thisDist - rad1 < 1.0) thisChar = 255 - (unsigned char)(255.0*0.5*(1.0 + thisDist - rad1));
        //printf("%d %d %d  %g %d\n",i,j,k,thisDist,thisChar);

      } else {
        // do it the new way, accounting for linearly-varying radius along the segment
        thisDist = min_dist_with_rad(x1,y1,z1,rad1, x2,y2,z2,rad2, (double)i+0.5,(double)j+0.5,(double)k+0.5);
      }

      // convert that distance to an unsigned char
      // 255=inside
      // 127=right on the boundary
      // 0=more than 1 voxel away from surface
      unsigned char thisChar = 0;
      if (thisDist < -1.0) thisChar = 255;
      else if (thisDist < 1.0) thisChar = 255 - (unsigned char)(255.0*0.5*(1.0 + thisDist));

      // only update the array if this voxel is nearer to this segment
      if (thisChar > dat[i][j][k]) dat[i][j][k] = thisChar;
    }
    }
    }

    // advance
    curr = curr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }
  fprintf(stderr,"\n");
  fflush(stderr);

  // finally, write the file
  fprintf(stderr,"  writing...\n");
  fflush(stderr);
  (void) write_bob_file_from_uchar(ofp, dat, nx, ny, nz);

  // free the memory and return
  free_3d_array_b(dat, nx, ny, nz);

  return(0);
}


/*
 * Read a 3D brick of bytes
 */
unsigned char*** readBobFile (char* inFileName, int* dims) {
  unsigned char*** data = NULL;
  int nx, ny, nz;
  FILE *fp;

  // open file for reading
  fp = fopen(inFileName,"rb");
  if (fp==NULL) {
    fprintf(stderr,"Could not open input file %s\n",inFileName);
    fprintf(stderr,"Exiting.\n");
    exit(1);
  }

  // read the header
  fread(&(nx),sizeof(int),1,fp);
  fread(&(ny),sizeof(int),1,fp);
  fread(&(nz),sizeof(int),1,fp);
  fprintf(stderr,"  reading bob file, size %d x %d x %d\n",nx,ny,nz);

  // make space for the array
  data = allocate_3d_array_b (nx, ny, nz);

  // read the file, plane-by-plane (safe)
  for (int i=0; i<nx; i++)
     fread(data[i][0],sizeof(char),ny*nz,fp);

  // assign grid dimensions
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;

  return data;
}

// needed externs
seg_ptr add_segment (seg_group_ptr, node_ptr, node_ptr);
node_ptr add_node (node_group_ptr, unsigned char, double*, int);

/*
 * Convert a 3D brick of bytes into segmented "3d graph paper" based on values
 */
int read_bob (char* infile, seg_group_ptr thisSG) {

  // read the data
  int dims[3] = {-1, -1, -1};
  unsigned char*** dat = readBobFile (infile, dims);
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];

  // recenter and rescale the data
  //int center = FALSE;
  float frac = 0.3;
  int scale = 14;
  int thresh = (int)(frac*256.0*scale*scale*scale);
  int nxnew = 1+nx/scale;
  int nynew = 1+ny/scale + 6;
  int nznew = 1+nz/scale;
  unsigned int*** scaled = allocate_3d_array_i (nxnew, nynew, nznew);
  for (int i=0; i<nxnew; i++)
  for (int j=0; j<nynew; j++)
  for (int k=0; k<nznew; k++) {
    scaled[i][j][k] = 0;
  }
  for (int i=0; i<nx; i++)
  for (int j=0; j<ny; j++)
  for (int k=0; k<nz; k++) {
    scaled[i/scale][j/scale][k/scale] += (int)dat[i][j][k];
  }
  free_3d_array_b(dat,nx,ny,nz);
  nx = nxnew;
  ny = nynew;
  nz = nznew;
  fprintf(stderr,"  bool array is %d x %d x %d\n",nx,ny,nz);

  // add a column
  // find the start
  if (FALSE) {
    int jlength = 6;
    int jstart = -1;
    int nxc = (nx-1)/2;
    for (int j=ny-1; j>-1; j--) {
      if (scaled[nxc][j][0] > thresh || scaled[nxc][j][1] > thresh) {
        jstart = j;
        break;
      }
    }
    if (jstart > -1) {
      for (int j=jstart; j<ny && j<jstart+jlength; j++) {
        scaled[nxc][j][0] = 2*thresh;
      }
    }
  }

  // booleanize the data
  nxnew = 1+nx;
  nynew = 1+ny;
  nznew = 1+nz;
  // "node on" and "cell on"
  unsigned char*** non = allocate_3d_array_b (nxnew, nynew, nznew);
  unsigned char*** con = allocate_3d_array_b (nxnew+1, nynew+1, nznew+1);
  for (int i=0; i<nxnew; i++) for (int j=0; j<nynew; j++) for (int k=0; k<nznew; k++) {
    non[i][j][k] = FALSE;
  }
  for (int i=0; i<nxnew+1; i++) for (int j=0; j<nynew+1; j++) for (int k=0; k<nznew+1; k++) {
    con[i][j][k] = FALSE;
  }
  int cellcnt = 0;
  for (int i=0; i<nx; i++)
  for (int j=0; j<ny; j++)
  for (int k=0; k<nz; k++) {
    if (scaled[i][j][k] > thresh) {
      cellcnt++;
      con[i+1][j+1][k+1] = TRUE;
      non[i][j][k] = TRUE;
      non[i][j][k+1] = TRUE;
      non[i][j+1][k] = TRUE;
      non[i][j+1][k+1] = TRUE;
      non[i+1][j][k] = TRUE;
      non[i+1][j][k+1] = TRUE;
      non[i+1][j+1][k] = TRUE;
      non[i+1][j+1][k+1] = TRUE;
    }
  }
  fprintf(stderr,"  found %d cells over thresh %d\n",cellcnt,thresh);
  free_3d_array_i(scaled,nx,ny,nz);
  nx = nxnew;
  ny = nynew;
  nz = nznew;

  // generate segments from boolean voxels
  thisSG->dim = 3;
  // box size is 1.0, default radius should be 1/(2r), where r is thickness ratio
  thisSG->radius = 0.07;
  node_ptr n1,n2;
  double tempd[3],otherd[3];
  for (int i=0; i<nx; i++) {
    tempd[0] = (double)i;
    for (int j=0; j<ny; j++) {
      tempd[1] = (double)j;
      for (int k=0; k<nz; k++) {
        tempd[2] = (double)k;
        if (non[i][j][k]) {
          //fprintf(stderr,"    node at %d %d %d\n",i,j,k);
          // make the node itself
          n1 = add_node (thisSG->nodes, thisSG->dim, tempd, 0);
          // connect segments
          if (i>0) if (non[i-1][j][k]) {
            if (con[i][j+1][k+1] || con[i][j+1][k] || con[i][j][k+1] || con[i][j][k]) {
              for (int ii=0; ii<3; ii++) otherd[ii] = tempd[ii];
              otherd[0] -= 1.0;
              n2 = add_node (thisSG->nodes, thisSG->dim, otherd, 0);
              (void) add_segment (thisSG, n1, n2);
            }
          }
          if (j>0) if (non[i][j-1][k]) {
            if (con[i+1][j][k+1] || con[i+1][j][k] || con[i][j][k+1] || con[i][j][k]) {
              for (int ii=0; ii<3; ii++) otherd[ii] = tempd[ii];
              otherd[1] -= 1.0;
              n2 = add_node (thisSG->nodes, thisSG->dim, otherd, 0);
              (void) add_segment (thisSG, n1, n2);
            }
          }
          if (k>0) if (non[i][j][k-1]) {
            if (con[i+1][j+1][k] || con[i+1][j][k] || con[i][j+1][k] || con[i][j][k]) {
              for (int ii=0; ii<3; ii++) otherd[ii] = tempd[ii];
              otherd[2] -= 1.0;
              n2 = add_node (thisSG->nodes, thisSG->dim, otherd, 0);
              (void) add_segment (thisSG, n1, n2);
            }
          }
        }
      }
    }
  }
  // wow, this is getting ugly, I think I need a proper 3d array object
  free_3d_array_b(non,nx,ny,nz);
  free_3d_array_b(con,nx+1,ny+1,nz+1);

  return (0);
}

