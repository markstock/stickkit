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
#include <png.h>
//#include "png.h"
#include "stickkit.h"

#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL

// file-local declarations
png_byte** allocate_2d_array_pb (int,int,int);
int free_2d_array_pb (png_byte**);
float** allocate_2d_array_f(int,int);
int free_2d_array_f (float**);

//void writePngPixels (const char*, const Array<png_byte,2>&, const bool);
//void readPngPixels (char*, png_byte**);
float** read_png_pixels (char*, float, float, int*, int*);

// externs from stickkit.c
extern node_ptr add_node (node_group_ptr, unsigned char, double*, int);
extern seg_ptr add_segment (seg_group_ptr, node_ptr, node_ptr);

/*
//
// scale and write an array to a png
//
int writeArrayToPng (Array<double,2>& data, int plotType, bool useBlackWhite, double black, double white) {

  static int frameCount = 0;

  // check the incoming Array
  if (data.dimensions() != 2) {
    cerr << "Error (writeArrayToPng): array not 2D!" << endl;
    exit(0);
  }

  // set the plot data range
  double blackval,whiteval;
  if (useBlackWhite) {
    blackval = black;
    whiteval = white;
  } else {
    blackval = min(data);
    whiteval = max(data);
  }

  // assemble file name
  char filename[16];
  sprintf(filename,"out_%05d.png",frameCount);
  cout << "    writing " << filename << " in range " << whiteval - blackval << endl << flush;

  // temporary array to store scaled data
  Array<double,2> scaled(data.shape());

  if (plotType == 3) {
    // illuminate only specific bands
    //float factor = 16.5 / (whiteval-blackval);
    //scaled = (data-blackval) * factor;
    scaled = linearToOneBand(data);
    scaled = scaled*255.0;

  } else if (plotType == 2) {
    // optionally scale to show even bands
    float factor = 16.5 / (whiteval-blackval);
    //float factor = 16.5 / 562.0;
    scaled = (data-blackval) * factor;
    scaled = linearToBands(scaled);
    scaled = scaled*255.0;

  } else if (plotType == 1) {
    // optionally scale to accentuate iso-lines
    float factor = 16.5 / (whiteval-blackval);
    scaled = (data-blackval) * factor;
    scaled = linearToIsolines(scaled);
    scaled = scaled*255.0;

  } else {
    // linearly scale the result
    float factor = 255.0 / (whiteval-blackval);
    scaled = (data-blackval) * factor;
  }

  // cast the array to the proper byte format
  Array<png_byte,2> bytes(data.shape());
  bytes = demoteToByte(scaled);
  scaled.free();

  // call the writer
  writePngPixels (filename,bytes,false);
  bytes.free();

  return frameCount++;
}


#ifndef __WIN32__
// write 8-bit color pixels to a PNG file
void writePngPixels (const char* filename, const Array<png_byte,2>& pxls, const bool isrgb) {

  int j;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  png_bytepp row_pointers = NULL;

  // we need to point into the Array data structure if we want to 
  //   avoid making a copy of the data
  row_pointers = (png_bytepp) malloc (pxls.rows() * sizeof (png_bytep));
  row_pointers[0] = (png_bytep)pxls.data();
  for (j=1; j<pxls.rows(); j++) {
    row_pointers[j] = row_pointers[j-1] + pxls.cols();
  }

  // open the file for writing
  FILE *fp = fopen(filename, "wb");
  if (fp==NULL) {
    fprintf(stderr,"Could not open output file (%s)\n",filename);
    fflush(stderr);
    exit(0);
  }

  // Create and initialize the png_struct with the desired error handler
  // functions.  If you want to use the default stderr and longjump method,
  // you can supply NULL for the last three parameters.  We also check that
  // the library version is compatible with the one used at compile time,
  // in case we are using dynamically linked libraries.  REQUIRED.
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (png_ptr == NULL) {
    fclose(fp);
    fprintf(stderr,"Could not create png struct\n");
    fflush(stderr);
    exit(0);
    // silent fail
    return;
  }

  // Allocate/initialize the image information data.  REQUIRED 
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fclose(fp);
    png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
    // silent fail
    return;
  }

  // Set error handling.  REQUIRED if you aren't supplying your own
  // error handling functions in the png_create_write_struct() call.
  if (setjmp(png_jmpbuf(png_ptr))) {
    // If we get here, we had a problem reading the file
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    // silent fail
    return;
  }

  // set up the output control if you are using standard C streams
  png_init_io(png_ptr, fp);

  // Set the image information here.  Width and height are up to 2^31,
  // bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
  // the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
  // PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
  // or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
  // PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
  // currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
  if (isrgb) {
  png_set_IHDR(png_ptr, info_ptr, (png_uint_32)pxls.cols(), (png_uint_32)pxls.rows(), 8,
    PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
    PNG_FILTER_TYPE_BASE);
  } else {
  png_set_IHDR(png_ptr, info_ptr, (png_uint_32)pxls.cols(), (png_uint_32)pxls.rows(), 8,
    PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
    PNG_FILTER_TYPE_BASE);
  }

  // Optional gamma chunk is strongly suggested if you have any guess
  // as to the correct gamma of the image.
  //png_set_gAMA(png_ptr, info_ptr, 2.2);
  //png_set_gAMA(png_ptr, info_ptr, 1.8);
  //png_set_gAMA(png_ptr, info_ptr, 0.55555);
  png_set_gAMA(png_ptr, info_ptr, 1.0);

  // Write the file header information.  REQUIRED
  png_write_info(png_ptr, info_ptr);

  // One of the following output methods is REQUIRED
  png_write_image(png_ptr, row_pointers);
  //png_write_image(png_ptr, imgrgb);

  // It is REQUIRED to call this to finish writing the rest of the file
  png_write_end(png_ptr, info_ptr);

  // clean up after the write, and free any memory allocated
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // close and increment
  fclose(fp);
}
#endif
*/


/*
 * read a PNG, convert it to segments along the given isoline
 */
int read_png(char* infile, seg_group_ptr this, double thresh) {

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
   int high_depth;
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
   if (bit_depth == 16) high_depth = TRUE;
   else high_depth = FALSE;

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
