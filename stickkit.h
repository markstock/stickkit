/*
 * stickkit.h
 *
 * copyright 2007,08,13 Mark J. Stock mstock@umich.edu
 *
 * a program to read, manipulate, and write string, segmented, or
 * graph-like data in fewer than 16 space dimensions
 *
 * compile with
 *
 * cc -O2 -std=c99 -o stickkit stickkit.c -lm
 * /usr/local/bin/i386-mingw32-gcc -O2 -o stickkit.exe stickkit.c -lm
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

// user-changeable parameters
#define BUCKET 100
#define DOTPER 10000
#define MAXDIM 16
#define MAX_ACTIONS 100

// necessary defines
#define TRUE 1
#define FALSE 0
#define EPSILON 1.0e-6
#define EPSILONSQRD 1.0e-12
#define M_PI           3.14159265358979323846
#define MAXSTR 160
#define DBLFLAG 1.234567d+99
//#define VERSION ("stickkit, version 0.1 2007-06-23 MJS")
//#define VERSION ("stickkit, version 0.2 2007-07-22 MJS")
//#define VERSION ("stickkit, version 0.3 2007-07-28 MJS")
// Version 0.4: added -info and -dots options
//#define VERSION ("stickkit, version 0.4 2008-04-29 MJS")
// Version 0.5: added -split option
//#define VERSION ("stickkit, version 0.5 2009-04-05 MJS")
// Version 0.6: added -rscale option
//#define VERSION ("stickkit, version 0.6 2014-03-06 MJS")
// Version 0.7: added -prune option
#define VERSION ("stickkit, version 0.7 2014-03-30 MJS")

// begin with the data structures

// predefine some structure pointers
typedef struct node *node_ptr;
typedef struct node_group *node_group_ptr;
typedef struct segment *seg_ptr;
typedef struct seg_group *seg_group_ptr;
typedef struct tangent *tan_ptr;
typedef struct tan_group *tan_group_ptr;
typedef struct radius *rad_ptr;
typedef struct rad_group *rad_group_ptr;

// a node
typedef struct node {
  unsigned long int index;
  char flag;
  double *x;
  unsigned int numconn0;
  seg_ptr *conn0;
  unsigned int numconn1;
  seg_ptr *conn1;
  node_group_ptr parent;
  node_ptr prev;
  node_ptr next;
} NODE;

// a collection of nodes, like a cell in a tree
typedef struct node_group {
  unsigned int num;
  node_ptr first;
  node_group_ptr parent;
  node_group_ptr child[2];
  unsigned char axis;
  double *min;
  double *max;
} NODE_GROUP;

// a segment
typedef struct segment {
  unsigned long int index;	// id of individual segment
  unsigned int block;		// id of disconnected block, start at 1
  char flag;			// temporary flag
  node_ptr n[2];
  rad_ptr r[2];
  tan_ptr t[2];
  seg_group_ptr parent;
  seg_ptr prev;
  seg_ptr next;
} SEGMENT;

// a collection of segments
typedef struct seg_group {
  unsigned char dim;		// dimensionality of system
  double radius;		// global radius, -1 if none defined
  unsigned long int num;	// number of segments in this group
  unsigned int numblock;	// number of segment blocks in this group
  seg_ptr first;		// pointer to first segment in linked list
  node_group_ptr nodes;		// pointer to group of nodes
  rad_group_ptr radii;		// pointer to group of radius objects
  tan_group_ptr tangents;	// pointer to group of tangent vectors
} SEGMENT_GROUP;

// a tangent vector
typedef struct tangent {
  unsigned long int index;
  char flag;
  double *x;
  tan_group_ptr parent;
  tan_ptr prev;
  tan_ptr next;
} TANGENT;

// a collection of tangent vectors
typedef struct tan_group {
  unsigned int num;
  tan_ptr first;
  tan_group_ptr parent;
  tan_group_ptr child[2];
  unsigned char axis;
  double *min;
  double *max;
} TANGENT_GROUP;

// a radius
typedef struct radius {
  unsigned long int index;
  char flag;
  double r;
  rad_group_ptr parent;
  rad_ptr prev;
  rad_ptr next;
} RADIUS;

// a collection of radii
typedef struct rad_group {
  unsigned int num;
  rad_ptr first;
  rad_group_ptr parent;
  rad_group_ptr child[2];
  double min;
  double max;
} RADIUS_GROUP;

// define each action possible on the data
typedef enum action_type {
  none,		// no value, no action (default)
  sk_read,	// read a file
  sk_write,	// write a file
  info,		// dump min/max bounds, number of nodes, segments
  coarsen,	// remove nodes and merge segments
  refine,	// split segments and add nodes, assume straight-line segments
  splrefine,	// split segments and add nodes, use spline interpolation
  roughen,	// nudge nodes to roughen the curve
  treeradius,	// set radii according to principle of constant stress
  translate,	// translate all nodes by x,y,z
  scale,	// scale all nodes by x,y,z
  rscale,	// scale all radii by either value 1, or linearly between values 1 and 2
  split,	// split all segments along the longest axis, write both to .rad files
  prune		// clip segments within a distance from a tip
  //resample,	// resample a strand to a fixed segment length
  //smooth,	// nudge nodes to smooth the curve
} ACTION_NAME;

typedef struct an_action {
  ACTION_NAME type;
  unsigned char nc;
  char **carg;
  unsigned char ni;
  int *iarg;
  unsigned char nd;
  double *darg;
} ACTION;

// define the possible output file types
typedef enum output_format_type {
   noout,	// default is no output
   seg,		// Wavefront .obj-like .seg file
   rad,		// radiance cylinders
   vtk,		// Mayavi-readable
   dots,	// ASCII space-delimited, equidistant dots, for use with laser subsurface etching
   png,		// image of segments, for debugging only
   svg,		// vector-based file, using wxSVG
   bob		// rasterized brick of bytes
} OUT_FORMAT;

