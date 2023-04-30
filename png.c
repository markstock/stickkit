/*
 * stickkit.c
 * png.c
 *
 * copyright 2013 Mark J. Stock mstock@umich.edu
 *
 * a program to read, manipulate, and write string, segmented, or
 * graph-like data in fewer than 16 space dimensions
 *
 * compile with
 *
 * g++ -O2 -c svgout.cpp
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
#include <math.h>
#include <png.h>
#include "stickkit.h"

#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL

// file-local declarations
png_byte** allocate_2d_array_pb (int,int,int);
int free_2d_array_pb (png_byte**);
float** allocate_2d_array_f(int,int);
int free_2d_array_f (float**);

float** read_png_pixels (char*, float, float, int*, int*);

int write_png_file (FILE*, const int, const int, const bool, float**, float, float);
int write_png (FILE*, seg_group_ptr, int);

// externs from stickkit.c
extern node_ptr add_node (node_group_ptr, unsigned char, double*, int);
extern seg_ptr add_segment (seg_group_ptr, node_ptr, node_ptr);

// externs from bob.c
extern int min(int, int);
extern int max(int, int);

// ------------- OUTPUT ----------------------------------------------------

//
// scale and write an array to a png
//
int write_png_file (FILE* fp, const int nx, const int ny, const bool high_depth,
   float **red, float redmin, float redrange) {

   bool autorange = true;
   int i,j,printval,bit_depth;
   char outfile[80];
   float newminrange,newmaxrange;
   // gamma of 1.8 looks normal on most monitors...that display properly.
   //float gamma = 1.8;
   // must do 5/9 for stuff to look right on Macs....why? I dunno.
   float gamma = .55555;
   //FILE *fp;
   png_uint_32 height,width;
   png_structp png_ptr;
   png_infop info_ptr;
   static png_byte **img;
   static bool is_allocated = false;

   // set specific bit depth
   if (high_depth) bit_depth = 16;
   else bit_depth = 8;

   // allocate the space for the special array
   if (!is_allocated) {
      img = allocate_2d_array_pb(nx,ny,bit_depth);
      is_allocated = true;
   }

   // set the sizes in png-understandable format
   height=ny;
   width=nx;

   // auto-set the ranges
   if (autorange) {

      // first red
      newminrange = 9.9e+9;
      newmaxrange = -9.9e+9;
      for (i=0; i<nx; i++) {
         for (j=ny-1; j>=0; j--) {
            if (red[i][j]<newminrange) newminrange=red[i][j];
            if (red[i][j]>newmaxrange) newmaxrange=red[i][j];
         }
      }
      //printf("range %g %g\n",newminrange,newmaxrange);
      redmin = newminrange;
      redrange = newmaxrange-newminrange;

   } else {
       // report the range
      newminrange = 9.9e+9;
      newmaxrange = -9.9e+9;
      for (i=0; i<nx; i++) {
         for (j=ny-1; j>=0; j--) {
            if (red[i][j]<newminrange) newminrange=red[i][j];
            if (red[i][j]>newmaxrange) newmaxrange=red[i][j];
         }
      }
      printf("  output range %g %g\n",newminrange,newmaxrange);
   }
 
   // make the preliminary filename (write the pgm to outfile)
   //sprintf(outfile,"%s.png",outfileroot);

   // write the file
   //fp = fopen(outfile,"wb");
   //if (fp==NULL) {
   //   fprintf(stderr,"Could not open output file %s\n",outfile);
   //   fflush(stderr);
   //   exit(0);
   //}

   // monochrome image, read data from red array
   // no scaling, 16-bit per channel
   if (high_depth) {
     // am I looping these coordinates in the right memory order?
     for (j=ny-1; j>=0; j--) {
       for (i=0; i<nx; i++) {
         printval = (int)(0.5 + 65535*(red[i][j]-redmin)/redrange);
         if (printval<0) printval = 0;
         else if (printval>65535) printval = 65535;
         img[ny-1-j][2*i] = (png_byte)(printval/256);
         img[ny-1-j][2*i+1] = (png_byte)(printval%256);
       }
     }

   // no scaling, 8-bit per channel
   } else {
     // am I looping these coordinates in the right memory order?
     for (j=ny-1; j>=0; j--) {
       for (i=0; i<nx; i++) {
         printval = (int)(0.5 + 256*(red[i][j]-redmin)/redrange);
         if (printval<0) printval = 0;
         else if (printval>255) printval = 255;
         img[ny-1-j][i] = (png_byte)printval;
       }
     }
   }

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);

   if (png_ptr == NULL) {
      //fclose(fp);
      fprintf(stderr,"Could not create png struct\n");
      fflush(stderr);
      exit(0);
      return (-1);
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      //fclose(fp);
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      return (-1);
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr))) {
      /* If we get here, we had a problem reading the file */
      //fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return (-1);
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth,
      PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
      PNG_FILTER_TYPE_BASE);

   /* Optional gamma chunk is strongly suggested if you have any guess
    * as to the correct gamma of the image. */
   //png_set_gAMA(png_ptr, info_ptr, 2.2);
   png_set_gAMA(png_ptr, info_ptr, gamma);

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* One of the following output methods is REQUIRED */
   png_write_image(png_ptr, img);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   // close file
   //fclose(fp);

   // free the data array
   free_2d_array_pb(img);

   return(0);
}

double min_dist_with_rad_2d (double vx, double vy, double vr,
                          double wx, double wy, double wr,
                          double px, double py) {
  // Return minimum distance between fat line segment vw and point p
  // i.e. |w-v|^2 -  avoid a sqrt
  const double l2 = pow(vx-wx,2) + pow(vy-wy,2);
  const double dl = sqrt(l2);
  const double dvp = sqrt( pow(vx-px,2) + pow(vy-py,2) );
  const double dwp = sqrt( pow(wx-px,2) + pow(wy-py,2) );
  // new: one radius overwhelms the entire segment and other endcap
  if (vr > dl+wr) return dvp - vr;
  if (wr > dl+vr) return dwp - wr;
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  double t = ( (px-vx)*(wx-vx) + (py-vy)*(wy-vy) ) / l2;
  double tadj = t;
  if (fabs(vr-wr) / dl > 1.e-5) {
    // segment is enough of a cone to matter
    // must shift t to account for cone
    // project to centerline, regardless of t
    const double jx = vx + t * (wx - vx);
    const double jy = vy + t * (wy - vy);
    const double pt = sqrt( pow(jx-px,2) + pow(jy-py,2) );
    // push t one direction or another
    tadj = t + (pt/dl) * (wr-vr) / sqrt(l2 - pow(wr-vr,2));
  }
  // Beyond the 'v' end of the segment
  if (tadj < 0.0) return dvp - vr;
  //if (tadj < 0.0) return dvp;
  // Beyond the 'w' end of the segment
  else if (tadj > 1.0) return dwp - wr;
  //else if (tadj > 1.0) return dwp;
  // Projection falls on the segment
  const double jx = vx + tadj * (wx - vx);
  const double jy = vy + tadj * (wy - vy);
  // and find the local radius
  const double thisRad = vr + tadj*(wr-vr);
  return sqrt( pow(jx-px,2) + pow(jy-py,2) ) - thisRad;
}

int write_png (FILE* ofp, seg_group_ptr thisSG, int res) {

  int xdim = 0;
  //int ydim = thisSG->dim-1;
  int ydim = 1;
  int nx,ny;
  double start[2];
  double size[2];

  // wrong if only 2D
  if (thisSG->dim > 2) {
    fprintf(stderr,"Will only print first and last dimensions to png.\n");
    fflush(stderr);
  } else if (thisSG->dim == 1) {
    fprintf(stderr,"How did you get a 1D .seg file?\n");
    fflush(stderr);
    return(1);
  }

  // eventually, find a smart way to project DIM dimensions down to 2

  // how big is this domain?
  //(void) update_node_group_stats(thisSG->nodes, thisSG->dim);
  fprintf(stderr,"  node minima");
  for (int i=0; i<thisSG->dim; i++) fprintf(stderr," %g",thisSG->nodes->min[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  node maxima");
  for (int i=0; i<thisSG->dim; i++) fprintf(stderr," %g",thisSG->nodes->max[i]);
  fprintf(stderr,"\n");

  // set build domain 10% larger (will this be enough to capture wide lines?)
  //for (int i=0; i<2; i++) size[i] = thisSG->nodes->max[i] - thisSG->nodes->min[i];
  size[0] = thisSG->nodes->max[xdim] - thisSG->nodes->min[xdim];
  size[1] = thisSG->nodes->max[ydim] - thisSG->nodes->min[ydim];
  double maxRad = thisSG->radii->max;
  if (thisSG->radii->max < 0.0) maxRad = thisSG->radius;
  fprintf(stderr,"  radius maximum %g\n",maxRad);
  //for (int i=0; i<2; i++) start[i] = thisSG->nodes->min[i] - maxRad;
  start[0] = thisSG->nodes->min[xdim] - maxRad;
  start[1] = thisSG->nodes->min[ydim] - maxRad;
  for (int i=0; i<2; i++) size[i] += 2.0*maxRad;

  // compute reasonable bounds
  double dx;
  if (size[0] > size[1]) {
    nx = res;
    dx = size[0] / nx;
    ny = 1 + (int)(size[1] / dx);
  } else {
    ny = res;
    dx = size[1] / ny;
    nx = 1 + (int)(size[0] / dx);
  }
  fprintf(stderr,"  cell size %g\n",dx);
  fprintf(stderr,"  png will be %d x %d\n",nx,ny);

  // sanity check on bob size
  if (nx*(float)ny > 1.e+10 || nx > 100000 || ny > 100000) {
    fprintf(stderr,"Are you sure you want to write a png file that large?\n");
    fprintf(stderr,"  Hit enter to continue, control-c to quit.\n");
    fflush(stderr);
    char buf[2];
    fread(&buf, sizeof(char), 1, stdin);
  }

  // allocate space for the brick
  float** dat = allocate_2d_array_f(nx,ny);

  // zero out
  for (int i=0; i<nx; i++) {
  for (int j=0; j<ny; j++) {
    dat[i][j] = 0;
  }
  }

  fprintf(stderr,"  computing");
  fflush(stderr);

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
    const double x1 = (curr->n[0]->x[xdim] - start[0]) / dx;
    const double y1 = (curr->n[0]->x[ydim] - start[1]) / dx;
    const double x2 = (curr->n[1]->x[xdim] - start[0]) / dx;
    const double y2 = (curr->n[1]->x[ydim] - start[1]) / dx;

    // find x,y,z range affected by this segment
    const int imin = max((int)floor(fmin(x1-rad1, x2-rad2)) - 1, 0);
    const int imax = min((int)ceil(fmax(x1+rad1, x2+rad2)) + 1, nx);
    const int jmin = max((int)floor(fmin(y1-rad1, y2-rad2)) - 1, 0);
    const int jmax = min((int)ceil(fmax(y1+rad1, y2+rad2)) + 1, ny);
    //printf("kmax  %d %d %d\n",(int)ceil(fmax(z1+rad1, z2+rad2)), nz, kmax);

    // loop over that subblock
    for (int i=imin; i<imax; i++) {
    for (int j=jmin; j<jmax; j++) {

      // how far is this node from the segment, in voxels?
      double thisDist = 2.0;

      thisDist = min_dist_with_rad_2d (x1,y1,rad1, x2,y2,rad2, (double)i+0.5,(double)j+0.5);

      // convert that distance to an unsigned char
      // 255=inside
      // 127=right on the boundary
      // 0=more than 1 voxel away from surface
      float thisVal = 0.0;
      if (thisDist < -1.0) thisVal = 1.0;
      else if (thisDist < 1.0) thisVal = 0.5 - 0.5*thisDist;

      // only update the array if this voxel is nearer to this segment
      //if (thisVal > dat[i][j]) dat[i][j] = thisVal;
      // or, always accumulate
      dat[i][j] += thisVal;
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

  (void) write_png_file (ofp, nx, ny, true, dat, 0.0, -1.0);

  // free and clear
  (void) free_2d_array_f (dat);

  return(0);
}


// -------------- INPUT ----------------------------------------------------

/*
 * read a PNG, convert it to segments along the given isoline
 */
int read_png (char* infile, seg_group_ptr this, double thresh) {

  // Image data
  int imgx, imgy;

  // Segment data
  double x[MAXDIM];
  int n0index,n1index;
  node_ptr newnode0 = NULL;
  node_ptr newnode1 = NULL;
  node_ptr *thenodes;
  seg_ptr newseg;

  fprintf(stderr,"  reading %s\n",infile);
  fflush(stderr);

  // Get the data from the image file
  float** dat = read_png_pixels(infile, 0.0, 1.0, &imgx, &imgy);

  // Prepare the data structure
  this->dim = 2;              // provide a default
  this->radius = 0.25;         // negative == "unset"
  this->num = 0;              // gotta start somewhere
  int nnode = 0;
  int nlines = 0;

  // Quickly count the number of nodes
  for (int j=0; j<imgy; j++) {
    for (int i=0; i<imgx-1; i++) {
      if ( (dat[i][j]-thresh) * (dat[i+1][j]-thresh) < 0.0) {
        nnode++;
      }
    }
  }
  for (int j=0; j<imgy-1; j++) {
    for (int i=0; i<imgx; i++) {
      if ( (dat[i][j]-thresh) * (dat[i][j+1]-thresh) < 0.0) {
        nnode++;
      }
    }
  }
  // Now, the number of lines
  for (int j=0; j<imgy-1; j++) {
    for (int i=0; i<imgx-1; i++) {
      int numCrossings = 0;

      // does the left edge contain a threshold-crossing?
      if ( (dat[i][j]-thresh) * (dat[i][j+1]-thresh) < 0.0) {
        numCrossings++;
      }

      // does the bottom edge contain a threshold-crossing?
      if ( (dat[i][j]-thresh) * (dat[i+1][j]-thresh) < 0.0) {
        numCrossings++;
      }

      // does the right edge contain a threshold-crossing?
      if ( (dat[i+1][j]-thresh) * (dat[i+1][j+1]-thresh) < 0.0) {
        numCrossings++;
      }

      // does the top edge contain a threshold-crossing?
      if ( (dat[i+1][j+1]-thresh) * (dat[i][j+1]-thresh) < 0.0) {
        numCrossings++;
      }

      nlines += numCrossings/2;
    }
  }
  fprintf(stderr,"  to read %d nodes, %d segments\n",nnode,nlines);
  fflush(stderr);
  nnode = 0;
  nlines = 0;

  // malloc space for an array of pointers
  //thenodes = (node_ptr*) malloc ((nnode+1) * sizeof(node_ptr));

  // Iterate over the pixels and generate segments, 0-2 per pixel
  for (int j=0; j<imgy-1; j++) {
    for (int i=0; i<imgx-1; i++) {

      // does the left edge contain a threshold-crossing?
      if ( (dat[i][j]-thresh) * (dat[i][j+1]-thresh) < 0.0) {
        x[0] = (double)i;
        x[1] = imgy - (double)j - (thresh-dat[i][j]) / (dat[i][j+1]-dat[i][j]);
        if (newnode0) {
          newnode1 = add_node (this->nodes, this->dim, x, 0);
          // create the segment
          newseg = add_segment (this, newnode0, newnode1);
          newnode0 = NULL;
        } else {
          newnode0 = add_node (this->nodes, this->dim, x, 0);
        }
      }

      // does the bottom edge contain a threshold-crossing?
      if ( (dat[i][j]-thresh) * (dat[i+1][j]-thresh) < 0.0) {
        x[0] = (double)i + (thresh-dat[i][j]) / (dat[i+1][j]-dat[i][j]);
        x[1] = imgy - (double)j;
        if (newnode0) {
          newnode1 = add_node (this->nodes, this->dim, x, 0);
          newseg = add_segment (this, newnode0, newnode1);
          newnode0 = NULL;
        } else {
          newnode0 = add_node (this->nodes, this->dim, x, 0);
        }
      }

      // does the right edge contain a threshold-crossing?
      if ( (dat[i+1][j]-thresh) * (dat[i+1][j+1]-thresh) < 0.0) {
        x[0] = (double)(i+1);
        x[1] = imgy - (double)j - (thresh-dat[i+1][j]) / (dat[i+1][j+1]-dat[i+1][j]);
        if (newnode0) {
          newnode1 = add_node (this->nodes, this->dim, x, 0);
          newseg = add_segment (this, newnode0, newnode1);
          newnode0 = NULL;
        } else {
          newnode0 = add_node (this->nodes, this->dim, x, 0);
        }
      }

      // does the top edge contain a threshold-crossing?
      if ( (dat[i+1][j+1]-thresh) * (dat[i][j+1]-thresh) < 0.0) {
        x[0] = (double)i + (thresh-dat[i][j+1]) / (dat[i+1][j+1]-dat[i][j+1]);
        x[1] = imgy - (double)(j+1);
        if (newnode0) {
          newnode1 = add_node (this->nodes, this->dim, x, 0);
          newseg = add_segment (this, newnode0, newnode1);
          newnode0 = NULL;
        } else {
          newnode0 = add_node (this->nodes, this->dim, x, 0);
        }
      }

    }
  }
  //fprintf(stderr,"%s has %d segments\n",infile,nlines);
  //fflush(stderr);

  return nlines;
}


/*
 * read a PNG, scale it, save it to 1 channel
 *
 * dat, nx, and ny are outputs
 */
float** read_png_pixels (char *infile, float datmin, float datrange, int* nx, int* ny) {

   int i,j;
   bool high_depth = false;
   FILE *fp;
   unsigned char header[8];
   float **dat;
   png_uint_32 height,width;
   int bit_depth,color_type,interlace_type;
   png_structp png_ptr;
   png_infop info_ptr;
   png_byte **img;

   // check the file
   fp = fopen(infile,"rb");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      fflush(stderr);
      exit(0);
   }

   // check to see that it's a PNG
   fread (&header, 1, 8, fp);
   if (png_sig_cmp(header, 0, 8)) {
      fprintf(stderr,"File %s is not a PNG\n",infile);
      fflush(stderr);
      exit(0);
   }

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also supply the
    * the compiler header file version, so that we know if the application
    * was compiled with a compatible version of the library.  REQUIRED
    */
   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);

   /* Allocate/initialize the memory for image information.  REQUIRED. */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      fclose(fp);
      png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
      exit(0);
   }

   /* Set error handling if you are using the setjmp/longjmp method (this is
    * the normal method of doing things with libpng).  REQUIRED unless you
    * set up your own error handlers in the png_create_read_struct() earlier.  */
   if (setjmp(png_jmpbuf(png_ptr))) {
      /* Free all of the memory associated with the png_ptr and info_ptr */
      png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
      fclose(fp);
      /* If we get here, we had a problem reading the file */
      exit(0);
   }

   /* One of the following I/O initialization methods is REQUIRED */
   /* Set up the input control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* If we have already read some of the signature */
   png_set_sig_bytes(png_ptr, 8);

   /* The call to png_read_info() gives us all of the information from the
    * PNG file before the first IDAT (image data chunk).  REQUIRED */
   png_read_info(png_ptr, info_ptr);

   png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
       &interlace_type, int_p_NULL, int_p_NULL);

   /* Set up the data transformations you want.  Note that these are all
    * optional.  Only call them if you want/need them.  Many of the
    * transformations only work on specific types of images, and many
    * are mutually exclusive.  */

   /* tell libpng to strip 16 bit/color files down to 8 bits/color */
   //png_set_strip_16(png_ptr);

   /* Extract multiple pixels with bit depths of 1, 2, and 4 from a single
    * byte into separate bytes (useful for paletted and grayscale images).  */
   png_set_packing(png_ptr);

   /* Expand paletted colors into true RGB triplets */
   //if (color_type == PNG_COLOR_TYPE_PALETTE)
   //   png_set_palette_rgb(png_ptr);

   /* Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel */
   if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
      png_set_expand_gray_1_2_4_to_8(png_ptr);

   /* Optional call to gamma correct and add the background to the palette
    * and update info structure.  REQUIRED if you are expecting libpng to
    * update the palette for you (ie you selected such a transform above).
    */
   //png_read_update_info(png_ptr, info_ptr);

   // check image type for applicability
   if (bit_depth != 8 && bit_depth != 16) {
     fprintf(stderr,"INCOMPLETE: read_png expect 8-bit or 16-bit images\n");
     fprintf(stderr,"   bit_depth: %d\n",bit_depth);
     fprintf(stderr,"   file: %s\n",infile);
     exit(0);
   }
   if (color_type != PNG_COLOR_TYPE_GRAY && color_type != PNG_COLOR_TYPE_RGB) {
     fprintf(stderr,"INCOMPLETE: read_png expect grayscale or RGB images\n");
     fprintf(stderr,"   color_type: %d\n",color_type);
     fprintf(stderr,"   file: %s\n",infile);
     exit(0);
   }

   // set channels
   if (color_type != PNG_COLOR_TYPE_GRAY) {
     fprintf(stderr,"ERROR: not expecting 3-channel PNG, but input is 3-channel\n");
     fprintf(stderr,"  file (%s)",infile);
     fprintf(stderr,"  Convert file to grayscale and try again.\n");
     exit(0);
   }

   // set specific bit depth
   if (bit_depth == 16) high_depth = true;
   else high_depth = false;

   // allocate the space for the image array
   img = allocate_2d_array_pb(width,height,bit_depth);

   /* Now it's time to read the image.  One of these methods is REQUIRED */
   png_read_image(png_ptr, img);

   /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
   png_read_end(png_ptr, info_ptr);

   /* At this point you have read the entire image */

   /* clean up after the read, and free any memory allocated - REQUIRED */
   png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

   /* close the file */
   fclose(fp);

   // allocate space for the 2D array of floats
   dat = allocate_2d_array_f(width,height);

   // monochrome image, read data from red array

   // no scaling, 16-bit per channel
   if (high_depth) {
      for (j=height-1; j>=0; j--) {
         for (i=0; i<width; i++) {
            dat[i][j] = datmin+datrange*(img[height-1-j][2*i]*256+img[height-1-j][2*i+1])/65535.;
         }
      }

   // no scaling, 8-bit per channel
   } else {
      for (j=height-1; j>=0; j--) {
         for (i=0; i<width; i++) {
            dat[i][j] = datmin+datrange*img[height-1-j][i]/255.;
         }
      }
   }

   // free the data array
   free_2d_array_pb(img);

   // set the sizes so that we can understand them
   (*nx) = width;
   (*ny) = height;

   return(dat);
}

/*
 * allocate memory for a two-dimensional array of png_byte
 */
png_byte** allocate_2d_array_pb(int nx, int ny, int depth) {

   int i,bytesperpixel;
   png_byte **array;

   if (depth <= 8) bytesperpixel = 1;
   else bytesperpixel = 2;
   array = (png_byte **)malloc(ny * sizeof(png_byte *));
   array[0] = (png_byte *)malloc(bytesperpixel * nx * ny * sizeof(png_byte));

   for (i=1; i<ny; i++)
      array[i] = array[0] + i * bytesperpixel * nx;

   return(array);
}

int free_2d_array_pb(png_byte** array){
   free(array[0]);
   free(array);
   return(0);
}

/*
 * allocate memory for a two-dimensional array of floats
 */
float** allocate_2d_array_f(int nx,int ny) {

   int i;
   float **array = (float **)malloc(nx * sizeof(float *));

   array[0] = (float *)malloc(nx * ny * sizeof(float));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}

int free_2d_array_f(float** array){
   free(array[0]);
   free(array);
   return(0);
}
