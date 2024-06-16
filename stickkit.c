/*
 * stickkit.c
 *
 * copyright 2007,08,13-15,17,24 Mark J. Stock mstock@umich.edu
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

// necessary defines
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "stickkit.h"
//#include <malloc.h>

// prototypes
// processing
void allocate_action (ACTION*, ACTION_NAME);
seg_group_ptr make_seg_structure (int);
int update_all_stats (seg_group_ptr);
int update_node_group_stats (node_group_ptr, int);
int update_radius_group_stats (rad_group_ptr);
int update_tangent_group_stats (tan_group_ptr, int);
int dump_stats (seg_group_ptr);
int prune_once (seg_group_ptr);
int merge_close_nodes (node_group_ptr, int, double);
void delete_segment (node_ptr, node_ptr, seg_ptr, int);
int split_long_segs (seg_group_ptr, double, int, int);
void split_segment (seg_ptr, int, int);
void find_seg_midpt_using_midpoint (double*, double*, seg_ptr, int);
void find_seg_midpt_using_spline  (double*, double*, seg_ptr, int);
void find_seg_tangents (seg_ptr, int);
int set_seg_tangents (seg_group_ptr, int);
int roughen_nodes (node_group_ptr, int, double);
int identify_separate_strands (seg_group_ptr);
int identify_connected_segments (seg_ptr, unsigned int);
int find_root_of_each_strand (seg_group_ptr, ACTION*);
//int set_treelike_radii_2d (seg_group_ptr, ACTION*);
int set_treelike_radii_3d (seg_group_ptr, ACTION*);
int translate_nodes (node_group_ptr, int, double*);
int scale_nodes (node_group_ptr, int, double*);
int scale_radii (rad_group_ptr, double*);
int split_and_write_rad (seg_group_ptr, char*);
int split_into_two_seg_groups (node_group_ptr, int, double, seg_group_ptr, seg_group_ptr);
// input/output
int read_seg (char*, seg_group_ptr, int);
int read_obj (char*, seg_group_ptr, int);
int read_rad (char*, seg_group_ptr);
int read_radiance_floats (FILE*, double*);
seg_ptr add_segment (seg_group_ptr, node_ptr, node_ptr);
node_ptr add_node (node_group_ptr, unsigned char, double*, int);
rad_ptr add_radius (rad_group_ptr, double, int);
tan_ptr add_tangent (tan_group_ptr, unsigned char, double*, int);
node_ptr find_closest_node (node_group_ptr, unsigned char, double *,
  node_ptr, double *);
int write_seg (FILE*, seg_group_ptr, int, char**);
int write_vtk (FILE*, seg_group_ptr);
int write_vtk_nodes_and_set_indexes (FILE*, node_group_ptr, int, int);
int write_vtk_nodes_radii (FILE*, node_group_ptr, double);
int write_svg (FILE*, seg_group_ptr, int, char**);
int write_rad (FILE*, seg_group_ptr);
int write_radiance_nodes (FILE*, node_group_ptr, int, double, int*, int*);
int write_dots (FILE*, seg_group_ptr, double);
int write_dots_1 (node_group_ptr, node_group_ptr, int, double *);
int write_dots_2 (seg_group_ptr, node_group_ptr, int, double *);
int write_dots_3 (FILE *, node_group_ptr, int);
int write_obj (FILE*, seg_group_ptr, double);
// utility
int set_seg_flags (seg_group_ptr, const bool);
int set_node_flags (node_group_ptr, const bool);
int set_radius_flags (rad_group_ptr, const bool);
void perturb_gaussian (double*, double, int);
double seg_length (seg_ptr, int);
double seg_length_sqrd (seg_ptr, int);
double vec_dist_sqrd (double*, double*, int);
double vec_inv_sqrt (double*, int);
double vec_dot (double const * const, double const * const, const int);
void vec_make_perpendicular (double * const, double const * const, const int);
void vec_normalize (double*, double*,const int);
void vec_normalize_in_place ( double*, const int);
void vec_cross (double*, double*, double*);
int free_nodes (node_group_ptr);
int free_radii (rad_group_ptr);
int free_tangents (tan_group_ptr);
int free_segs (seg_group_ptr);
int Usage(char[MAXSTR],int);

// in other files
extern int write_svg_using_wxSVG(FILE*, seg_group_ptr, int, char**);
//extern int writeArrayToPng (Array<double,2>&,int,bool,double,double);
//extern seg_group_ptr readArrayFromPng (char*, double**);
extern int read_png(char*, seg_group_ptr, double);
extern int write_bob(FILE*, seg_group_ptr, double);
extern int write_png(FILE*, seg_group_ptr, int);
extern int read_bob (char*, seg_group_ptr);

// do all the things
int main (int argc, char **argv) {

  // main non-global definitions
  int i,j,act;
  bool have_input_file = false;
  double dot_scale = 0.3;
  double bob_cell = 1.0;
  double obj_min = 1.0;
  int png_res = 1024;
  bool zero_indexed = false;
  seg_group_ptr segs;
  //seg_group_ptr segs = (SEGMENT_GROUP*) malloc (sizeof(SEGMENT_GROUP));
  char progname[MAXSTR],infile[MAXSTR],root[MAXSTR],extension[4];
  OUT_FORMAT out_format = noout;
  int nactions = 0;
  ACTION action[MAX_ACTIONS];
  FILE *outptr = stdout;

  // parse command-line arguments
  // input file can appear anywhere
  // processing options are performed in the order they appear
  (void) strcpy(progname,argv[0]);
  if (argc < 2) (void) Usage(progname,0);
  for (i=1; i<argc; i++) {
    if (strncmp(argv[i], "-",1) == 0) {
      if (strncmp(argv[i], "-seg", 4) == 0) {
        out_format = seg;
      } else if (strncmp(argv[i], "-svg", 4) == 0) {
        out_format = svg;
      } else if (strncmp(argv[i], "-rad", 4) == 0) {
        out_format = rad;
      } else if (strncmp(argv[i], "-vtk", 4) == 0) {
        out_format = vtk;
      } else if (strncmp(argv[i], "-dots", 4) == 0) {
        out_format = dots;
        if (argc > i+1) {
          if (!isalpha(argv[i+1][1])) {
            dot_scale = atof(argv[++i]);
          }
        }
      } else if (strncmp(argv[i], "-png", 4) == 0) {
        out_format = png;
        if (argc > i+1) {
          if (!isalpha(argv[i+1][1])) {
            png_res = atoi(argv[++i]);
          }
        }
      } else if (strncmp(argv[i], "-bob", 4) == 0) {
        out_format = bob;
        if (argc > i+1) {
          if (!isalpha(argv[i+1][1])) {
            bob_cell = atof(argv[++i]);
          }
        }
      } else if (strncmp(argv[i], "-obj", 4) == 0) {
        out_format = obj;
        if (argc > i+1) {
          if (!isalpha(argv[i+1][1])) {
            obj_min = atof(argv[++i]);
          }
        }
      } else if (strncmp(argv[i], "-noop", 5) == 0) {
        (void) allocate_action (&action[nactions], none);
        nactions++;
      } else if (strncmp(argv[i], "-coarsen", 5) == 0) {
        (void) allocate_action (&action[nactions], coarsen);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
        }
        nactions++;
      } else if (strncmp(argv[i], "-srefine", 5) == 0) {
        (void) allocate_action (&action[nactions], splrefine);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
        }
        nactions++;
      } else if (strncmp(argv[i], "-refine", 5) == 0) {
        (void) allocate_action (&action[nactions], refine);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
        }
        nactions++;
      } else if (strncmp(argv[i], "-roughen", 5) == 0) {
        (void) allocate_action (&action[nactions], roughen);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
        }
        nactions++;
      } else if (strncmp(argv[i], "-prune", 3) == 0) {
        (void) allocate_action (&action[nactions], prune);
        if (action[nactions].iarg[0] = atoi(argv[++i])) {
          // if we could convert it to an int check its value
          if (action[nactions].iarg[0] < 1) {
            action[nactions].iarg[0] = 1;
          }
          nactions++;
        } else {
          // if it wasn't an int, or it was 0, just don't prune
        }
      } else if (strncmp(argv[i], "-gr", 3) == 0) {
        (void) allocate_action (&action[nactions], globalradius);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
        }
        nactions++;
      } else if (strncmp(argv[i], "-treeradius", 4) == 0) {
        (void) allocate_action (&action[nactions], treeradius);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 1.e+3;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
          if (isalpha(argv[i+1][1])) {
            action[nactions].iarg[0] = 0;
          } else {
            action[nactions].iarg[0] = atoi(argv[++i]);
            if (isalpha(argv[i+1][1])) {
              action[nactions].darg[1] = 0.;
            } else {
              action[nactions].darg[1] = atof(argv[++i]);
            }
          }
        }
        nactions++;
      } else if (strncmp(argv[i], "-translate", 4) == 0) {
        (void) allocate_action (&action[nactions], translate);
        if (isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = 0.;
        } else {
          action[nactions].darg[0] = atof(argv[++i]);
          if (isalpha(argv[i+1][1])) {
            action[nactions].darg[1] = 0.;
          } else {
            action[nactions].darg[1] = atof(argv[++i]);
            if (isalpha(argv[i+1][1])) {
              action[nactions].darg[2] = 0.;
            } else {
              action[nactions].darg[2] = atof(argv[++i]);
            }
          }
        }
        nactions++;
      } else if (strncmp(argv[i], "-scale", 3) == 0) {
        (void) allocate_action (&action[nactions], scale);
        for (j=0;j<MAXDIM;j++) action[nactions].darg[j] = 1.;
        if (argc > i+1) {
          if (!isalpha(argv[i+1][1])) {
            i++;
            for (j=0;j<MAXDIM;j++) action[nactions].darg[j] = atof(argv[i]);
            if (argc > i+1) {
              if (!isalpha(argv[i+1][1])) {
                i++;
                for (j=1;j<MAXDIM;j++) action[nactions].darg[j] = atof(argv[i]);
                if (argc > i+1) {
                  if (!isalpha(argv[i+1][1])) {
                    i++;
                    for (j=2;j<MAXDIM;j++) action[nactions].darg[j] = atof(argv[i]);
                  }
                }
              }
            }
          }
        }
        nactions++;
      } else if (strncmp(argv[i], "-rscale", 3) == 0) {
        (void) allocate_action (&action[nactions], rscale);
        action[nactions].darg[0] = DBLFLAG;
        action[nactions].darg[1] = DBLFLAG;
        if (!isalpha(argv[i+1][1])) {
          action[nactions].darg[0] = atof(argv[++i]);
          if (!isalpha(argv[i+1][1])) {
            action[nactions].darg[1] = atof(argv[++i]);
          }
        }
        nactions++;
      } else if (strncmp(argv[i], "-split", 3) == 0) {
        (void) allocate_action (&action[nactions], split);
        nactions++;
      } else if (strncmp(argv[i], "-zeroindex", 2) == 0) {
        zero_indexed = true;
      } else if (strncmp(argv[i], "-info", 2) == 0) {
        (void) allocate_action (&action[nactions], info);
        nactions++;
      } else if (strncmp(argv[i], "--V", 3) == 0) {
        fprintf(stderr,"You ran (%s)\n",progname);
        fprintf(stderr,"%s\n",VERSION);
        exit(0);
      } else if (strncmp(argv[i], "-version", 2) == 0) {
        fprintf(stderr,"You ran (%s)\n",progname);
        fprintf(stderr,"%s\n",VERSION);
        exit(0);
      } else {
        fprintf(stderr,"Unrecognized option (%s)\n\n",argv[i]);
        (void) Usage(progname,0);
      }
    } else {
      // arg does not start with "-" must be an input filename!
      //(void) strcpy (infile,argv[i]);
      have_input_file = true;
      (void) allocate_action (&action[nactions], sk_read);
      (void) strcpy (action[nactions].carg[0], argv[i]);
      nactions++;
    }
  }

  // set things up, if necessary

  // last action is a write if a file format has been specified
  if (out_format != none) {
    (void) allocate_action (&action[nactions], sk_write);
    nactions++;
  }

  // eventally support reading in from stdin!
  if (!have_input_file) {
    fprintf(stderr,"ERROR: no input file present on command-line.\n");
    fprintf(stderr,"Exiting\n");
    exit(0);
  }

  // prepare the segs structure, argument is dimensions
  segs = make_seg_structure (3);

  // now, march through the list of actions, performing one at a time
  for (act=0; act<nactions; act++) {

    fprintf (stderr,"Action %d of %d : ",act+1,nactions);

    // read one or more files from the file lists
    if (action[act].type == sk_read) {

      fprintf (stderr,"read\n");
      // check extension
      (void) strcpy (infile,action[act].carg[0]);
      //fprintf(stderr,"(%s)\n",infile); fflush(stderr);
      strncpy (extension,infile+strlen(infile)-3,4);
      // read the input file
      if (strncmp(extension, "seg", 3) == 0) {
        (void) read_seg (infile, segs, zero_indexed);
      } else if (strncmp(extension, "obj", 3) == 0) {
        (void) read_obj (infile, segs, zero_indexed);
      } else if (strncmp(extension, "rad", 3) == 0) {
        (void) read_rad (infile, segs);
      } else if (strncmp(extension, "bob", 3) == 0) {
        (void) read_bob (infile, segs);
      } else if (strncmp(extension, "png", 3) == 0) {
        (void) read_png (infile, segs, 0.5);
      } else {
        fprintf (stderr,"ERROR: input file format (%s) is not",extension);
        fprintf (stderr," supported yet.\nSupported file types: seg,obj,rad,png\n");
        fprintf (stderr,"Exiting\n");
        exit(0);
      }

    // write the new strands
    } else if (action[act].type == sk_write) {

      fprintf (stderr,"write\n");
      if (out_format == svg) {
        (void) write_svg (outptr, segs, argc, argv);
      } else if (out_format == vtk) {
        (void) write_vtk (outptr, segs);
      } else if (out_format == rad) {
        (void) write_rad (outptr, segs);
      } else if (out_format == dots) {
        (void) write_dots (outptr, segs, dot_scale);
      } else if (out_format == obj) {
        (void) write_obj (outptr, segs, obj_min);
      } else if (out_format == bob) {
        (void) write_bob (outptr, segs, bob_cell);
      } else if (out_format == png) {
        (void) write_png (outptr, segs, png_res);
      } else {
        (void) write_seg (outptr, segs, argc, argv);
      }

    // remove nodes closer than a threshold distance
    } else if (action[act].type == coarsen) {

      fprintf (stderr,"coarsen\n");
      fprintf (stderr,"  begin with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);
      // do the merge, interate until all are done
      i = 1;
      while (i > 0) {
        (void) set_node_flags (segs->nodes, false);
        i = merge_close_nodes (segs->nodes, segs->dim,
                                  pow(action[act].darg[0],2));
        fprintf (stderr,"  performed %d merges\n",i);
      }
      // reset the counters and bounds
      (void) update_all_stats (segs);
      fprintf (stderr,"  end with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);

    // remove nodes closer than a threshold distance
    } else if (action[act].type == splrefine) {

      fprintf (stderr,"spline refine\n");
      fprintf (stderr,"  begin with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);

      // should we do the normalized version? if thresh is negative
      if (action[act].darg[0] > 0.) j = false;
      else j = true;
      //fprintf (stderr,"  normalization is %d\n",j);

      // do the splits, interate until all are done
      i = 1;
      while (i > 0) {
        (void) set_seg_flags (segs, false);
        i = split_long_segs (segs, pow(action[act].darg[0],2), true, j);
        fprintf (stderr,"  performed %d splits\n",i);
      }
      // reset the counters and bounds
      (void) update_all_stats (segs);
      fprintf (stderr,"  end with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);

    // remove nodes closer than a threshold distance
    } else if (action[act].type == refine) {

      fprintf (stderr,"refine\n");
      fprintf (stderr,"  begin with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);

      // should we do the normalized version? if thresh is negative
      if (action[act].darg[0] > 0.) j = false;
      else j = true;
      //fprintf (stderr,"  normalization is %d\n",j);

      // do the splits, interate until all are done
      i = 1;
      while (i > 0) {
        (void) set_seg_flags (segs, false);
        i = split_long_segs (segs, pow(action[act].darg[0],2), false, j);
        fprintf (stderr,"  performed %d splits\n",i);
      }
      // reset the counters and bounds
      (void) update_all_stats (segs);
      fprintf (stderr,"  end with %d segs, %d nodes\n",
               segs->num,segs->nodes->num);

    // jiggle nodes around
    } else if (action[act].type == roughen) {

      fprintf (stderr,"roughen\n");
      (void) roughen_nodes (segs->nodes, segs->dim, action[act].darg[0]);

    // set the global radius
    } else if (action[act].type == globalradius) {

      fprintf (stderr,"global radius set to %g\n", action[act].darg[0]);
      segs->radius = action[act].darg[0];

    // set radii according to principle of constant stress
    } else if (action[act].type == treeradius) {

      fprintf (stderr,"treeradius\n");
      if (segs->dim == 3) {
        // first, identify separate disconnected strands
        i = identify_separate_strands (segs);
        fprintf (stderr,"  found %d separate strands\n",i);
        // then, find the root node of each strand
        i = find_root_of_each_strand (segs, &action[act]);
        fprintf (stderr,"  found %d strand roots\n",i);
        // finally, set the radii
        (void) set_seg_flags (segs, false);
        i = 1;
        while (i > 0) {
          //i = set_treelike_radii_2d (segs, &action[act]);
          i = set_treelike_radii_3d (segs, &action[act]);
          fprintf (stderr,".");
        }
        fprintf (stderr,"\n");

      } else {
        fprintf (stderr,"  Error: treeradius cannot work on systems with");
        fprintf (stderr," %d dimensions\n",segs->dim);
        fprintf (stderr,"  works on dim 2 and 3\n");
        fprintf (stderr,"  skipping this step\n");
      }

    // remove nodes closer than a threshold distance
    } else if (action[act].type == translate) {

      fprintf (stderr,"translate\n");
      (void) translate_nodes (segs->nodes, segs->dim, action[act].darg);

    // scale all nodes
    } else if (action[act].type == scale) {

      fprintf (stderr,"scale\n");
      (void) scale_nodes (segs->nodes, segs->dim, action[act].darg);

    // scale all radii
    } else if (action[act].type == rscale) {

      fprintf (stderr,"scale radii\n");
      (void) scale_radii (segs->radii, action[act].darg);

    // split and write .rad files
    } else if (action[act].type == split) {

      // find input filename root
      if (strrchr(infile,'.')) {
        strncpy(root,infile,strrchr(infile,'.')-infile);
        root[strrchr(infile,'.')-infile] = '_';
        root[strrchr(infile,'.')-infile+1] = '\0';
      } else {
        // there is no "." in the filename?
        sprintf(root,"%s_",infile);
      }
      //fprintf(stderr,"root is (%s)\n",root);

      fprintf (stderr,"split\n");
      (void) split_and_write_rad (segs, root);

    // remove nodes and segments close to tips
    } else if (action[act].type == prune) {

      fprintf (stderr,"prune\n");
      for (i=0; i<action[act].iarg[0]; i++) {
        int num_pruned = prune_once (segs);
        fprintf (stderr,"  pruned %d segments\n",num_pruned);
      }

    // dump statistics and quit
    } else if (action[act].type == info) {

      fprintf (stderr,"info\n");
      (void) dump_stats (segs);

    // no action
    } else {

      fprintf(stderr,"none\n");
    }
  }

  exit(0);
}


/*
 * malloc the necessary space for this type of action
 */
void allocate_action (ACTION *action, ACTION_NAME type) {

  int i;

  action->type = type;

  if (type == sk_read) {
    // a file read, allocate one character variable, none others
    action->nc = 1;
    action->ni = 0;
    action->nd = 0;
  } else if (type == sk_write) {
    // a file write, allocate one character variable, none others
    action->nc = 1;
    action->ni = 0;
    action->nd = 0;
  } else if (type == info) {
    // dump current segments information
    action->nc = 0;
    action->ni = 0;
    action->nd = 0;
  } else if (type == coarsen) {
    // need a threshold
    action->nc = 0;
    action->ni = 0;
    action->nd = 1;
  } else if (type == splrefine) {
    // need a threshold
    action->nc = 0;
    action->ni = 0;
    action->nd = 1;
  } else if (type == refine) {
    // need a threshold
    action->nc = 0;
    action->ni = 0;
    action->nd = 1;
  } else if (type == roughen) {
    // need a scale
    action->nc = 0;
    action->ni = 0;
    action->nd = 1;
  } else if (type == treeradius) {
    // need a vector
    action->nc = 0;
    action->ni = 1;
    action->nd = 2;
  } else if (type == translate) {
    // need a vector
    action->nc = 0;
    action->ni = 0;
    action->nd = MAXDIM;
  } else if (type == scale) {
    // need a vector
    action->nc = 0;
    action->ni = 0;
    action->nd = MAXDIM;
  } else if (type == rscale) {
    // need a vector
    action->nc = 0;
    action->ni = 0;
    action->nd = 2;
  } else if (type == split) {
    // no need for arguments
    action->nc = 0;
    action->ni = 0;
    action->nd = 0;
  } else if (type == prune) {
    // no need for arguments
    action->nc = 0;
    action->ni = 1;
    action->nd = 0;
  } else if (type == globalradius) {
    // need a vector
    action->nc = 0;
    action->ni = 0;
    action->nd = 1;
  } else if (type == none) {
    // no operation!
    action->nc = 0;
    action->ni = 0;
    action->nd = 0;
  } else {
    fprintf (stderr,"ERROR (allocate_action): not supposed to get here.\n");
    fprintf (stderr,"Did you forget to edit allocate_action for the new operation?.\n");
    fprintf (stderr,"Quitting.\n");
    exit(0);
  }

  // now, malloc the space
  if (action->nc > 0) {
    action->carg = (char**) malloc (action->nc * sizeof(char*));
    for (i=0; i<action->nc; i++)
      action->carg[i] = (char*) malloc (MAXSTR * sizeof(char));
  } else {
    action->carg = NULL;
  }
  if (action->ni > 0) {
    action->iarg = (int*) malloc (action->ni * sizeof(int));
  } else {
    action->iarg = NULL;
  }
  if (action->nd > 0) {
    action->darg = (double*) malloc (action->nd * sizeof(double));
  } else {
    action->darg = NULL;
  }

  return;
}


/*
 * Create an independent set of related nodes and segments
 *
 * Should we support having multiple segment groups active?
 */
seg_group_ptr make_seg_structure (int dim) {

  seg_group_ptr segs = (SEGMENT_GROUP*) malloc (sizeof(SEGMENT_GROUP));

  segs->dim = dim;
  segs->num = 0;
  segs->first = NULL;
  segs->nodes = (NODE_GROUP*) malloc (sizeof(NODE_GROUP));
  segs->nodes->num = 0;
  segs->nodes->first = NULL;
  segs->nodes->parent = NULL;
  segs->nodes->child[0] = NULL;
  segs->nodes->child[1] = NULL;
  segs->nodes->axis = -1;
  segs->nodes->min = NULL;
  segs->nodes->max = NULL;
  segs->radii = (RADIUS_GROUP*) malloc (sizeof(RADIUS_GROUP));
  segs->radii->num = 0;
  segs->radii->first = NULL;
  segs->radii->parent = NULL;
  segs->radii->child[0] = NULL;
  segs->radii->child[1] = NULL;
  segs->radii->min = 9.9e+9;
  segs->radii->max = -9.9e+9;
  segs->tangents = (TANGENT_GROUP*) malloc (sizeof(TANGENT_GROUP));
  segs->tangents->num = 0;
  segs->tangents->first = NULL;
  segs->tangents->parent = NULL;
  segs->tangents->child[0] = NULL;
  segs->tangents->child[1] = NULL;
  segs->tangents->axis = -1;
  segs->tangents->min = NULL;
  segs->tangents->max = NULL;

  return (segs);
}


// ==========================================================================
// process the data
//

/*
 * Update numbers, bounds of all data
 */
int update_all_stats (seg_group_ptr thisSG) {

  int nrad, nnode, nseg; // ntan
  seg_ptr curr;

  // do segments first (they're easy)
  nseg = 0;
  curr = thisSG->first;
  while (curr) {
    curr->parent = thisSG;
    nseg++;
    curr = curr->next;
  }
  thisSG->num = nseg;

  // then do nodes
  nnode = update_node_group_stats (thisSG->nodes, thisSG->dim);

  // do radii next
  nrad = update_radius_group_stats (thisSG->radii);

  // do tangents last
  //ntan = update_tangent_group_stats (thisSG->tangents, thisSG->dim);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Print the system statistics
 */
int dump_stats (seg_group_ptr thisSG) {

  int i;

  fprintf(stdout,"# number of dimensions %d\n",thisSG->dim);
  fprintf(stdout,"# number of segment blocks %d\n",thisSG->numblock+1);
  fprintf(stdout,"# number of segments %ld\n",thisSG->num);
  fprintf(stdout,"# number of nodes %ld\n",thisSG->nodes->num);
  fprintf(stdout,"#   node minima");
  for (i=0; i<thisSG->dim; i++) fprintf(stdout," %g",thisSG->nodes->min[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"#   node maxima");
  for (i=0; i<thisSG->dim; i++) fprintf(stdout," %g",thisSG->nodes->max[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"# number of radii %ld\n",thisSG->radii->num);
  if (thisSG->radii->num > 0) {
    fprintf(stdout,"#   radius minimum %g\n",thisSG->radii->min);
    fprintf(stdout,"#   radius maximum %g\n",thisSG->radii->max);
  }
  fprintf(stdout,"# number of tangents %ld\n",thisSG->tangents->num);

  return(0);
}


/*
 * Update the num, min, max of a group of nodes
 */
int update_node_group_stats (node_group_ptr thisNG, int dim) {

  int i,j;
  int cnt = 0;
  node_ptr curr;

  // reset the min/max
  for (i=0; i<dim; i++) {
    thisNG->min[i] = 9.9e+9;
    thisNG->max[i] = -9.9e+9;
  }

  // check all children
  if (thisNG->child[0]) {

    // update counts and bounds of children
    cnt += update_node_group_stats (thisNG->child[0], dim);
    cnt += update_node_group_stats (thisNG->child[1], dim);

    // use those bounds to update this bounds
    for (j=0; j<2; j++) {
      for (i=0; i<dim; i++) {
        if (thisNG->child[j]->min[i] < thisNG->min[i])
          thisNG->min[i] = thisNG->child[j]->min[i];
        if (thisNG->child[j]->max[i] > thisNG->max[i])
          thisNG->max[i] = thisNG->child[j]->max[i];
      }
    }

  } else {

    // check all local nodes
    curr = thisNG->first;
    while (curr) {

      // update parent
      curr->parent = thisNG;

      // update count
      cnt++;

      // update bounds
      for (i=0; i<dim; i++) {
        if (curr->x[i] < thisNG->min[i])
          thisNG->min[i] = curr->x[i];
        if (curr->x[i] > thisNG->max[i])
          thisNG->max[i] = curr->x[i];
      }

      curr = curr->next;
    }
  }

  thisNG->num = cnt;

  return (cnt);
}


/*
 * Update the num, min, max of a group of radii
 */
int update_radius_group_stats (rad_group_ptr thisRG) {

  int cnt = 0;
  rad_ptr curr;

  // if this group is NULL, skip and return
  if (thisRG) {

    // reset the min/max
    thisRG->min = 9.9e+9;
    thisRG->max = -9.9e+9;

    // first, count the local entries and compute min/max
    curr = thisRG->first;
    while (curr) {

      // update parent
      curr->parent = thisRG;

      // update count
      cnt++;

      // update bounds
      if (curr->r < thisRG->min) thisRG->min = curr->r;
      if (curr->r > thisRG->max) thisRG->max = curr->r;
      curr = curr->next;
    }

    // then, add on all children
    if (thisRG->child[0]) {
      cnt += update_radius_group_stats (thisRG->child[0]);
      // and check the children's bounds
      if (thisRG->child[0]->min < thisRG->min) thisRG->min = thisRG->child[0]->min;
      if (thisRG->child[0]->max > thisRG->max) thisRG->max = thisRG->child[0]->max;
    }
    if (thisRG->child[1]) {
      cnt += update_radius_group_stats (thisRG->child[1]);
      if (thisRG->child[1]->min < thisRG->min) thisRG->min = thisRG->child[1]->min;
      if (thisRG->child[1]->max > thisRG->max) thisRG->max = thisRG->child[1]->max;
    }

    // set the new count
    thisRG->num = cnt;
  }

  return (cnt);
}


/*
 * Merge nodes that share a common segment if their distance is under
 * a threshold (the threshold is squared already!)
 *
 * cnt is the number of merges that took place
 */
int merge_close_nodes (node_group_ptr thisNG, int dim, double thresh) {

  int i;
  int cnt = 0;
  double dist,mindist;
  node_ptr curr,delnode;
  seg_ptr delseg;

  // check all children
  if (thisNG->child[0]) {
    cnt += merge_close_nodes (thisNG->child[0], dim, thresh);
    cnt += merge_close_nodes (thisNG->child[1], dim, thresh);

  } else {
    // check all connected segments
    //fprintf (stderr,"group with %d nodes, thresh %g\n",thisNG->num,thresh);

    curr = thisNG->first;
    while (curr) {

      // loop through all connected segments
      //fprintf (stderr, "node %ld at %g %g %g has neighbors at\n",
      //         curr->index, curr->x[0], curr->x[1], curr->x[2]);

      // loop until there are no more close nodes
      //mindist = 0.;
      //while (mindist < thresh) {

      // no, do looping higher up!

        // find the closest adjacent connected node
        mindist = 9.9e+9;
        delnode = NULL;
        delseg = NULL;

        for (i=0; i<curr->numconn0; i++) {
          // how close is the connected node?
          dist = vec_dist_sqrd (curr->x, curr->conn0[i]->n[1]->x, dim);

          //fprintf (stderr, "  %ld at %g %g %g dist %g\n",
          //         curr->conn0[i]->n[1]->index,curr->conn0[i]->n[1]->x[0],
          //         curr->conn0[i]->n[1]->x[1],curr->conn0[i]->n[1]->x[2],
          //         sqrt(dist));

          if (dist < mindist) {
            mindist = dist;
            delnode = curr->conn0[i]->n[1];
            delseg = curr->conn0[i];
          }
        }

        for (i=0; i<curr->numconn1; i++) {
          dist = vec_dist_sqrd (curr->x, curr->conn1[i]->n[0]->x, dim);

          //fprintf (stderr, "  %ld at %g %g %g dist %g\n",
          //         curr->conn1[i]->n[0]->index,curr->conn1[i]->n[0]->x[0],
          //         curr->conn1[i]->n[0]->x[1],curr->conn1[i]->n[0]->x[2],
          //         sqrt(dist));

          if (dist < mindist) {
            mindist = dist;
            delnode = curr->conn1[i]->n[0];
            delseg = curr->conn1[i];
          }
        }

        // delete them both
        if (mindist < thresh && !curr->flag && !delnode->flag) {
          delete_segment (curr, delnode, delseg, dim);
          cnt++;
        }
      //}

      curr = curr->next;
    }
  }

  return (cnt);
}


/*
 * General routine for deleting a segment (and joining two nodes)
 */
void delete_segment (node_ptr keepnode, node_ptr delnode, seg_ptr delseg,
                     int dim) {

  int i,num_link_keep,num_link_del,copyto;
  int keepis0 = -1;
  double weight,newloc[dim];
  seg_ptr *newconnlist;
  rad_ptr newrad,keeprad,delrad;

  //fprintf (stderr, "     should we delete?\n");

  // first, make sure that the two given nodes are actually joined
  //   by the segment in question!
  if (delseg->n[0] == keepnode) {
    if (delseg->n[1] == delnode) {
      keepis0 = true;
    } else {
      keepis0 = -1;	// bad!
    }
  } else if (delseg->n[1] == keepnode) {
    if (delseg->n[0] == delnode) {
      keepis0 = false;
    } else {
      keepis0 = -1;	// bad!
    }
  } else {
    keepis0 = -1;	// bad!
  }

  // if they're not right, dump and quit
  if (keepis0 == -1) {
    fprintf (stderr,"ERROR (delete_segment): segment %ld was given nodes\n",delseg->index);
    fprintf (stderr,"  (%ld) and (%ld), but that segment's end nodes are\n",keepnode->index,delnode->index);
    fprintf (stderr,"  actually (%ld) and (%ld). Quitting.\n",delseg->n[0]->index,delseg->n[1]->index); 
    exit(1);
  }
  //fprintf (stderr, "    deleting node %ld and seg %ld, keeping node %ld\n",delnode->index,delseg->index,keepnode->index);

  // throw away the index and the flag of the old node

  // and set the new flag to indicate that it's participated this merge
  keepnode->flag = true;

  // relocate keepnode, weight by number of connections
  num_link_keep = keepnode->numconn0 + keepnode->numconn1;
  num_link_del = delnode->numconn0 + delnode->numconn1;
  //fprintf (stderr, "    num_link_keep %d, num_link_del %d\n",num_link_keep,num_link_del);
  if (num_link_keep == 1) {
    weight = 0.0;
  } else if (num_link_del == 1) {
    weight = 1.0;
  } else {
    weight = (double)(num_link_keep-1)/(double)(num_link_keep+num_link_del-2);
  }

  // This is fixed now---any simply-supported segment is cleanly deleted
  for (i=0; i<dim; i++) {
    keepnode->x[i] = weight*keepnode->x[i] + (1.-weight)*delnode->x[i];
  }
  // now, find the new radius, if there are radii at all
  if (delseg->parent->radii->first) {
    if (keepis0) {
      keeprad = delseg->r[0];
      delrad = delseg->r[1];
      newrad = add_radius (delseg->parent->radii, weight*delseg->r[0]->r +
                           (1.-weight)*delseg->r[1]->r, 0);
    } else {
      keeprad = delseg->r[1];
      delrad = delseg->r[0];
      newrad = add_radius (delseg->parent->radii, weight*delseg->r[1]->r +
                           (1.-weight)*delseg->r[0]->r, 0);
    }
  } else {
    keeprad = NULL;
    newrad = NULL;
  }

  // do something with this!
  // assign it to all still-connected segments
  for (i=0; i<keepnode->numconn0; i++)
    if (keepnode->conn0[i]->r[0] == keeprad)
      keepnode->conn0[i]->r[0] = newrad;
  for (i=0; i<keepnode->numconn1; i++)
    if (keepnode->conn1[i]->r[1] == keeprad)
      keepnode->conn1[i]->r[1] = newrad;
  for (i=0; i<delnode->numconn0; i++)
    if (delnode->conn0[i]->r[0] == delrad)
      delnode->conn0[i]->r[0] = newrad;
  for (i=0; i<delnode->numconn1; i++)
    if (delnode->conn1[i]->r[1] == delrad)
      delnode->conn1[i]->r[1] = newrad;

  // copy all conn0 linkages from delnode to keepnode
  //fprintf (stderr, "    numconn0: %d + %d",keepnode->numconn0,delnode->numconn0);
  if (delnode->numconn0 > 0) {
    newconnlist = (seg_ptr*) malloc (
                  (keepnode->numconn0 + delnode->numconn0) * sizeof(seg_ptr));
    for (i=0; i<keepnode->numconn0; i++) newconnlist[i] = keepnode->conn0[i];
    for (i=0; i<delnode->numconn0; i++) {
      delnode->conn0[i]->n[0] = keepnode;
      // add these to keepnode->conn0
      newconnlist[i + keepnode->numconn0] = delnode->conn0[i];
    }
    free (keepnode->conn0);
    keepnode->conn0 = newconnlist;
    keepnode->numconn0 += delnode->numconn0;
    //fprintf (stderr, " = %d\n",keepnode->numconn0);
  }

  // copy all conn1 linkages from delnode to keepnode
  //fprintf (stderr, "    numconn1: %d + %d",keepnode->numconn1,delnode->numconn1);
  if (delnode->numconn1 > 0) {
    newconnlist = (seg_ptr*) malloc (
                  (keepnode->numconn1 + delnode->numconn1) * sizeof(seg_ptr));
    for (i=0; i<keepnode->numconn1; i++) newconnlist[i] = keepnode->conn1[i];
    for (i=0; i<delnode->numconn1; i++) {
      delnode->conn1[i]->n[1] = keepnode;
      // add these to keepnode->conn1
      newconnlist[i + keepnode->numconn1] = delnode->conn1[i];
    }
    free (keepnode->conn1);
    keepnode->conn1 = newconnlist;
    keepnode->numconn1 += delnode->numconn1;
    //fprintf (stderr, " = %d\n",keepnode->numconn1);
  }

  // if it's the first node, reset the first pointer
  if (delnode->parent->first == delnode)
    delnode->parent->first = delnode->next;
  delnode->parent->num--;

  // reset pointers
  if (delnode->prev) delnode->prev->next = delnode->next;
  if (delnode->next) delnode->next->prev = delnode->prev;

  // remove delnode
  //fprintf (stderr, "    deleting node %ld\n",delnode->index);
  free (delnode->conn0);
  free (delnode->conn1);
  free (delnode);

  // now, remove delseg

  // if it's the first segment, reset the first pointer
  if (delseg->parent->first == delseg)
    delseg->parent->first = delseg->next;
  delseg->parent->num--;

  // lose everything except the pointers
  if (delseg->prev) delseg->prev->next = delseg->next;
  if (delseg->next) delseg->next->prev = delseg->prev;

  // remove *both* entries to delseg in keepnode's conn array
  // locate delseg in the conn0 array and compress the array
  copyto = 0;
  for (i=0; i<keepnode->numconn0; i++) {
    keepnode->conn0[copyto] = keepnode->conn0[i];
    if (keepnode->conn0[i] != delseg) copyto++;
  }
  keepnode->numconn0 = copyto;
  copyto = 0;
  // compress the conn1 array instead
  for (i=0; i<keepnode->numconn1; i++) {
    keepnode->conn1[copyto] = keepnode->conn1[i];
    if (keepnode->conn1[i] != delseg) copyto++;
  }
  keepnode->numconn1 = copyto;
  // don't worry about reallocating memory here

  // debug print the keepnode's connectivity list
  //fprintf (stderr, "    numconn0 are\n");
  //for (i=0; i<keepnode->numconn0; i++) {
  //  fprintf (stderr, "      seg %ld, nodes %ld %ld\n",keepnode->conn0[i]->index,keepnode->conn0[i]->n[0]->index,keepnode->conn0[i]->n[1]->index);
  //}
  //fprintf (stderr, "    numconn1 are\n");
  //for (i=0; i<keepnode->numconn1; i++) {
  //  fprintf (stderr, "      seg %ld, nodes %ld %ld\n",keepnode->conn1[i]->index,keepnode->conn1[i]->n[0]->index,keepnode->conn1[i]->n[1]->index);
  //}

  // finally, axe the segment
  free (delseg);

  //fprintf (stderr, "    done deleting\n");

  return;
}


/*
 * Split nodes that share a common segment if their distance is over
 * a threshold (the threshold is squared already!)
 *
 * cnt is the number of splits that took place
 * thresh is the *square* of the input threshold length
 */
int split_long_segs (seg_group_ptr thisSG, double thresh, int use_splines, int use_normalized_lengths) {

  int cnt = 0;
  double dist,rad;
  seg_ptr curr;

  curr = thisSG->first;
  while (curr) {

    // how long is this segment
    dist = vec_dist_sqrd (curr->n[0]->x, curr->n[1]->x, thisSG->dim);

    // scale the length by the radius
    if (use_normalized_lengths) {
      rad = 0.5*(curr->r[0]->r + curr->r[1]->r);
      dist /= pow(rad,2);
      //fprintf(stderr,"dist %g  rad %g  thresh %g  doit %d\n",sqrt(dist),rad,sqrt(thresh),dist>thresh);
    }

    // split this segment
    if (dist > thresh && !curr->flag) {
      split_segment (curr, thisSG->dim, use_splines);
      cnt++;
    }

    curr = curr->next;
  }

  return (cnt);
}


/*
 * General routine for splitting a segment
 */
void split_segment (seg_ptr thisSP, int dim, int use_spline) {

  int i,copyto;
  double radval,newloc[dim];
  seg_ptr newseg = NULL;
  node_ptr newnode = NULL;
  rad_ptr newrad = NULL;

  fprintf (stderr, "splitting seg %ld, nodes %ld %ld\n",thisSP->index,thisSP->n[0]->index,thisSP->n[1]->index);

  // locate the new node

  if (use_spline) {
    // spline-fit
    find_seg_midpt_using_spline (newloc, &radval, thisSP, dim);

  } else {
    // naive method (midpoint)
    find_seg_midpt_using_midpoint (newloc, &radval, thisSP, dim);
  }

  // make the new node
  newnode = add_node (thisSP->parent->nodes, dim, newloc, 0);

  // make the segment
  newseg = add_segment (thisSP->parent, newnode, thisSP->n[1]);

  // now, set the new radius
  newrad = add_radius (thisSP->parent->radii, radval, 0);

  // assign a bunch of pointers!
  // first, the new node
  newnode->numconn0 = 1;
  newnode->conn0 = (seg_ptr*) malloc (newnode->numconn0 * sizeof(seg_ptr));
  newnode->conn0[0] = newseg;
  newnode->numconn1 = 1;
  newnode->conn1 = (seg_ptr*) malloc (newnode->numconn1 * sizeof(seg_ptr));
  newnode->conn1[0] = thisSP;

  // but now this segment is no longer connected to its n[1]
  copyto = 0;
  for (i=0; i<thisSP->n[1]->numconn1; i++) {
    thisSP->n[1]->conn1[copyto] = thisSP->n[1]->conn1[i];
    if (thisSP->n[1]->conn1[i] != thisSP) copyto++;
  }
  thisSP->n[1]->numconn1 = copyto;

  // then the new segment
  newseg->r[0] = newrad;
  newseg->r[1] = thisSP->r[1];
  newseg->t[1] = thisSP->t[1];

  // and the old segment
  thisSP->n[1] = newnode;
  thisSP->r[1] = newrad;
  thisSP->t[1] = NULL;

  // set the two seg flags to indicate that they've participated in a split
  thisSP->flag = true;
  newseg->flag = true;

  return;
}


/*
 * use a naive method to determine the new node location - the midpoint
 */
void find_seg_midpt_using_midpoint (double *loc, double *rad,
    seg_ptr thisSP, int dim) {

  int i;
  double r0,r1;

  // midpoint position
  for (i=0;i<dim;i++) loc[i] = 0.5*(thisSP->n[0]->x[i] + thisSP->n[1]->x[i]);

  // midpoint radius
  if (thisSP->r[0]) r0 = thisSP->r[0]->r;
  else r0 = thisSP->parent->radius;
  if (thisSP->r[1]) r1 = thisSP->r[1]->r;
  else r1 = thisSP->parent->radius;
  *rad = 0.5*r0 + 0.5*r1;

  return;
}


/*
 * Make the new location at the midpoint of a cubic spline threaded
 *   through the two endpoints
 */
void find_seg_midpt_using_spline (double *loc, double *rad,
    seg_ptr thisSP, int dim) {

  int i,j;
  //double dl[dim],len,fp[2][3][3],p1[dim],p2[dim];
  double a[4];
  double dl,tan0[dim],tan1[dim];
  double r0,r1,r2;
  //double difflen = 0.;
  //static double maxdiff = 0.;
  //node_ptr n1 = thisSP->n[0];
  //node_ptr n2 = thisSP->n[1];

  // test print
  //fprintf(stderr,"\nnode0 at %g %g %g\n",thisSP->n[0]->x[0],thisSP->n[0]->x[1],thisSP->n[0]->x[2]);
  //fprintf(stderr,"node1 at %g %g %g\n",thisSP->n[1]->x[0],thisSP->n[1]->x[1],thisSP->n[1]->x[2]);

  // is a tangent defined here?
  fprintf(stderr,"making tangents\n");
  if (!thisSP->t[0] || !thisSP->t[1]) find_seg_tangents (thisSP, dim);

  //fprintf(stderr,"  tangent0 %g %g %g\n",thisSP->t[0]->x[0],thisSP->t[0]->x[1],thisSP->t[0]->x[2]);
  //fprintf(stderr,"  tangent1 %g %g %g\n",thisSP->t[1]->x[0],thisSP->t[1]->x[1],thisSP->t[1]->x[2]);

  // compute the length of the edge
/*
  len = 0.;
  for (i=0; i<dim; i++) {
    dl[i] = n2->x[i] - n1->x[i];
    len += dl[i]*dl[i];
  }
  len = sqrt(len);

  // compute the tangential operator for each node (P = I - nn^T)
  for (i=0; i<dim; i++) {
    for (j=0; j<3; j++) {
      fp[0][i][j] = 0.;
      fp[1][i][j] = 0.;
    }
    fp[0][i][i] = 1.;
    fp[1][i][i] = 1.;
    for (j=0; j<3; j++) {
      fp[0][i][j] -= n1->norm[i]*n1->norm[j];
      fp[1][i][j] -= n2->norm[i]*n2->norm[j];
    }
  }
  // find the vector product of each of these with dl
  for (i=0; i<dim; i++) {
    p1[i] = 0.;
    p2[i] = 0.;
  }
  for (i=0; i<dim; i++) {
    for (j=0; j<3; j++) {
      p1[i] += dl[j]*fp[0][i][j];
      p2[i] += dl[j]*fp[1][i][j];
    }
  }
*/

  // it looks like tangent curves should be scaled by dl
  dl = sqrt(vec_dist_sqrd (thisSP->n[0]->x, thisSP->n[1]->x, dim));
  for (i=0; i<dim; i++) tan0[i] = dl * thisSP->t[0]->x[i];
  for (i=0; i<dim; i++) tan1[i] = dl * thisSP->t[1]->x[i];

  // compute a spline for each coordinate axis
  for (i=0; i<dim; i++) {

    // a[0] = f_0
    a[0] = thisSP->n[0]->x[i];

    // a[1] = f'_0
    a[1] = tan0[i];

    // a[2] = 3 (f_1 - f_0) / h^2 - (f'_1 + 2 f'_0) / h
    // but h=1 for us
    a[2] = 3. * (thisSP->n[1]->x[i] - thisSP->n[0]->x[i])
              - (tan1[i] + 2. * tan0[i]);

    // a[3] = 2 (f_0 - f_1) / h^3 + (f'_1 + f'_0) / h^2
    a[3] = 2. * (thisSP->n[0]->x[i] - thisSP->n[1]->x[i])
              + (tan0[i] + tan1[i]);

    // and the node that we want is 1/2 along the curve
    loc[i] = a[0] + 0.5*a[1] + 0.25*a[2] + 0.125*a[3];

    //fprintf(stdout,"  a %g %g %g %g\n",a[0],a[1],a[2],a[3]);
  }

  //fprintf(stderr,"  real node loc %g %g %g\n",loc[0],loc[1],loc[2]);
  //exit(0);

  // Now, do the radius (must find slope of radius?) -------------------

  // find the end radii of this segment
  if (thisSP->r[0]) r0 = thisSP->r[0]->r;
  else r0 = thisSP->parent->radius;
  if (thisSP->r[1]) r1 = thisSP->r[1]->r;
  else r1 = thisSP->parent->radius;

  // find the length of this segment
  dl = sqrt(vec_dist_sqrd (thisSP->n[0]->x, thisSP->n[1]->x, dim));

  //fprintf(stderr,"\ndoing radius between %g and %g, length %g\n",r0,r1,dl);
  if (r0 < -EPSILON || r1 < -EPSILON) exit(0);

  // find the radius of the next node past n[0]
  j = thisSP->n[0]->numconn0 + thisSP->n[0]->numconn1;
  if (j == 1) {
    // n0 is an end, what to do?
    tan0[0] = r1-r0;
  } else if (j==2) {
    // there's one node past n[0], find its radius
    for (i=0; i<thisSP->n[0]->numconn0; i++) {
      if (thisSP->n[0]->conn0[i] != thisSP) {
        // find the length of that segment
        tan0[1] = sqrt (vec_dist_sqrd (thisSP->n[0]->conn0[i]->n[0]->x,
                                       thisSP->n[0]->conn0[i]->n[1]->x, dim));
        // find the radius of that segment
        if (thisSP->n[0]->conn0[i]->r[1]) r2 = thisSP->n[0]->conn0[i]->r[1]->r;
        else r2 = thisSP->parent->radius;
        // find the slope of the tangent
        //tan0[0] = (r1-r2) / (dl+tan0[1]);
        tan0[0] = 0.5*(r0-r2)*dl/tan0[1] + 0.5*(r1-r0);
      }
    }
    for (i=0; i<thisSP->n[0]->numconn1; i++) {
      if (thisSP->n[0]->conn1[i] != thisSP) {
        // find the length of that segment
        tan0[1] = sqrt (vec_dist_sqrd (thisSP->n[0]->conn1[i]->n[0]->x,
                                       thisSP->n[0]->conn1[i]->n[1]->x, dim));
        // find the radius of that segment
        if (thisSP->n[0]->conn1[i]->r[0]) r2 = thisSP->n[0]->conn1[i]->r[0]->r;
        else r2 = thisSP->parent->radius;
        // find the slope of the tangent
        //tan0[0] = (r1-r2) / (dl+tan0[1]);
        tan0[0] = 0.5*(r0-r2)*dl/tan0[1] + 0.5*(r1-r0);
      }
    }
    //fprintf(stderr,"  other seg ends with rad %g and len %g\n",r2,tan0[1]);
  } else {
    // if there are 3 segments connected to n[0], then just assume linear
    tan0[0] = r1-r0;
  }

  // find the radius of the next node past n[1]
  j = thisSP->n[1]->numconn0 + thisSP->n[1]->numconn1;
  if (j == 1) {
    // n0 is an end, what to do?
    tan1[0] = r1-r0;
  } else if (j==2) {
    // there's one node past n[1], find its radius
    for (i=0; i<thisSP->n[1]->numconn0; i++) {
      if (thisSP->n[1]->conn0[i] != thisSP) {
        // find the length of that segment
        tan1[1] = sqrt (vec_dist_sqrd (thisSP->n[1]->conn0[i]->n[0]->x,
                                       thisSP->n[1]->conn0[i]->n[1]->x, dim));
        // find the radius of that segment
        if (thisSP->n[1]->conn0[i]->r[1]) r2 = thisSP->n[1]->conn0[i]->r[1]->r;
        else r2 = thisSP->parent->radius;
        // find the slope of the tangent
        //tan1[0] = (r2-r0) / (dl+tan1[1]);
        tan1[0] = 0.5*(r2-r1)*dl/tan1[1] + 0.5*(r1-r0);
      }
    }
    for (i=0; i<thisSP->n[1]->numconn1; i++) {
      if (thisSP->n[1]->conn1[i] != thisSP) {
        // find the length of that segment
        tan1[1] = sqrt (vec_dist_sqrd (thisSP->n[1]->conn1[i]->n[0]->x,
                                       thisSP->n[1]->conn1[i]->n[1]->x, dim));
        // find the radius of that segment
        if (thisSP->n[1]->conn1[i]->r[0]) r2 = thisSP->n[1]->conn1[i]->r[0]->r;
        else r2 = thisSP->parent->radius;
        // find the slope of the tangent
        //tan1[0] = (r2-r0) / (dl+tan1[1]);
        tan1[0] = 0.5*(r2-r1)*dl/tan1[1] + 0.5*(r1-r0);
      }
    }
    //fprintf(stderr,"  other seg ends with rad %g and len %g\n",r2,tan1[1]);
  } else {
    // if there are 3 segments connected to n[1], then just assume linear
    tan1[0] = r1-r0;
  }

  //fprintf(stderr,"  tangents %g and %g\n",tan0[0],tan1[0]);

  // do the spline interpolation
  a[0] = r0;
  a[1] = tan0[0];
  a[2] = 3. * (r1 - r0) - (tan1[0] + 2. * tan0[0]);
  a[3] = 2. * (r0 - r1) + (tan0[0] + tan1[0]);
  *rad = a[0] + 0.5*a[1] + 0.25*a[2] + 0.125*a[3];
  //fprintf(stderr,"  new radius %g\n",*rad);
  if (*rad < EPSILON) {
    //fprintf(stderr,"  new radius %g\n",*rad);
    *rad = 0.;
    //exit(0);
  }

  return;
}


/*
 * Find the tangent vectors for a segment
 */
void find_seg_tangents (seg_ptr thisSP, int dim) {

  int i,j;
  int do_later0 = false;
  int do_later1 = false;
  double dotp,thisrad,otherrad,radsum;
  double x[dim],dl[dim],dlnorm[dim],other[dim];
  node_ptr n0 = thisSP->n[0];
  node_ptr n1 = thisSP->n[1];

  // we may need dl
  for (i=0; i<dim; i++) dl[i] = n1->x[i] - n0->x[i];
  fprintf(stderr,"  dl %g %g\n",dl[0],dl[1]);
  (void) vec_normalize (dl, dlnorm, dim);
  fprintf(stderr,"  dlnorm %g %g\n",dlnorm[0],dlnorm[1]);

  // do we have a tangent vector already?
  if (!thisSP->t[0]) {
    fprintf(stderr,"  no t[0] %d %d\n", n0->numconn0, n0->numconn1);

    // create the tangent
    thisSP->t[0] = add_tangent (thisSP->parent->tangents, dim, dl, 0);

    // if this is an end, then mirror the other tangent
    if (n0->numconn0 + n0->numconn1 == 1) {
      do_later0 = true;

    // otherwise, use summations
    } else {

      // now, loop over connected segments, counting this one
      for (i=0; i<dim; i++) x[i] = 0.;
      radsum = 0.;
      //fprintf(stderr,"    conn %d %d\n",n0->numconn0,n0->numconn1);
      for (j=0; j<n0->numconn0; j++) {
        // sum the radius-weighted, normalized tangent vectors
        for (i=0; i<dim; i++) other[i] = n0->conn0[j]->n[1]->x[i] -
                                         n0->conn0[j]->n[0]->x[i];
        (void) vec_normalize (other, other, dim);
        fprintf(stderr,"      other %g %g\n",other[0],other[1]);
        if (n0->conn0[j]->r[0])
          otherrad = pow (n0->conn0[j]->r[0]->r, 2);
        else
          otherrad = pow (thisSP->parent->radius, 2);
        fprintf(stderr,"      otherrad %g\n",otherrad);
        for (i=0; i<dim; i++) x[i] += other[i] * otherrad;
        radsum += otherrad;
      }
      for (j=0; j<n0->numconn1; j++) {
        for (i=0; i<dim; i++) other[i] = n0->conn1[j]->n[1]->x[i] -
                                         n0->conn1[j]->n[0]->x[i];
        (void) vec_normalize (other, other, dim);
        fprintf(stderr,"      other %g %g\n",other[0],other[1]);
        if (n0->conn1[j]->r[1])
          otherrad = pow (n0->conn1[j]->r[1]->r, 2);
        else
          otherrad = pow (thisSP->parent->radius, 2);
        fprintf(stderr,"      otherrad %g\n",otherrad);
        for (i=0; i<dim; i++) x[i] += other[i] * otherrad;
        radsum += otherrad;
      }
      fprintf(stderr,"    dlsum %g %g %g\n",x[0],x[1],x[2]);
      fprintf(stderr,"    radsum %g\n",radsum);

      if (thisSP->r[0]) thisrad = pow (thisSP->r[0]->r, 2);
      else thisrad = pow (thisSP->parent->radius, 2);
      // how important is this segment compared to all other segs?
      thisrad = 2.*thisrad/radsum;

      // finally, if this node just has 2 segments, use the same
      //   tangent on the other end
      if (n0->numconn0 + n0->numconn1 == 2) {
        fprintf(stderr,"    two-sided\n");

        // just use the weighted norm
        for (i=0; i<dim; i++) thisSP->t[0]->x[i] = x[i];
        (void) vec_normalize (thisSP->t[0]->x, thisSP->t[0]->x, dim);

        // and copy it to the other segment
        for (j=0; j<n0->numconn0; j++) {
          if (n0->conn0[j] != thisSP) {
            if (!n0->conn0[j]->t[0]) {
              fprintf(stderr,"      add_tangent0 %g %g\n",thisSP->t[0]->x[0],thisSP->t[0]->x[1]);
              n0->conn0[j]->t[0] = add_tangent (thisSP->parent->tangents,
                dim, thisSP->t[0]->x, 0);
        } } }
        for (j=0; j<n0->numconn1; j++) {
          if (n0->conn1[j] != thisSP) {
            if (!n0->conn1[j]->t[1]) {
              fprintf(stderr,"      add_tangent1 %g %g\n",thisSP->t[0]->x[0],thisSP->t[0]->x[1]);
              n0->conn1[j]->t[1] = add_tangent (thisSP->parent->tangents,
                dim, thisSP->t[0]->x, 0);
        } } }

      } else {

        // if this segment is tiny, then just plug it in straight
        // if its half, then just take the mean
        // if its a lot more, then, um, 
        for (i=0; i<dim; i++)
          thisSP->t[0]->x[i] = thisrad*x[i] + (1.-thisrad)*dlnorm[i];
        (void) vec_normalize (thisSP->t[0]->x, thisSP->t[0]->x, dim);
      }
    }

  }

  // is there one on this side?
  if (!thisSP->t[1]) {
    fprintf(stderr,"  no t[1] %d %d\n", n1->numconn0, n1->numconn1);

    thisSP->t[1] = add_tangent (thisSP->parent->tangents, dim, dl, 0);

    // if this is an end, then mirror the other tangent
    if (n1->numconn0 + n1->numconn1 == 1) {
      do_later1 = true;

    // otherwise, use summations
    } else {

      // now, loop over connected segments, counting this one
      for (i=0; i<dim; i++) x[i] = 0.;
      radsum = 0.;
      for (j=0; j<n1->numconn0; j++) {
        // sum the radius-weighted, normalized tangent vectors
        for (i=0; i<dim; i++) other[i] = n1->conn0[j]->n[1]->x[i] -
                                         n1->conn0[j]->n[0]->x[i];
        (void) vec_normalize (other, other, dim);
        if (n1->conn0[j]->r[0])
          otherrad = pow (n1->conn0[j]->r[0]->r, 2);
        else
          otherrad = pow (thisSP->parent->radius, 2);
        for (i=0; i<dim; i++) x[i] += other[i] * otherrad;
        radsum += otherrad;
      }
      for (j=0; j<n1->numconn1; j++) {
        for (i=0; i<dim; i++) other[i] = n1->conn1[j]->n[1]->x[i] -
                                         n1->conn1[j]->n[0]->x[i];
        (void) vec_normalize (other, other, dim);
        if (n1->conn1[j]->r[1])
          otherrad = pow (n1->conn1[j]->r[1]->r, 2);
        else
          otherrad = pow (thisSP->parent->radius, 2);
        for (i=0; i<dim; i++) x[i] += other[i] * otherrad;
        radsum += otherrad;
      }
      fprintf(stderr,"    dlsum %g %g %g\n",x[0],x[1],x[2]);
      fprintf(stderr,"    radsum %g\n",radsum);

      if (thisSP->r[1]) thisrad = pow (thisSP->r[1]->r, 2);
      else thisrad = pow (thisSP->parent->radius, 2);
      // how important is this segment compared to all other segs?
      thisrad = 2.*thisrad/radsum;

      // finally, if this node just has 2 segments, use the same
      //   tangent on the other end
      if (n1->numconn0 + n1->numconn1 == 2) {
        fprintf(stderr,"    two-sided\n");

        // just use the weighted norm
        for (i=0; i<dim; i++) thisSP->t[1]->x[i] = x[i];
        (void) vec_normalize (thisSP->t[1]->x, thisSP->t[1]->x, dim);

        // and copy it to the other segment
        for (j=0; j<n1->numconn0; j++) {
          if (n1->conn0[j] != thisSP) {
            if (!n1->conn0[j]->t[0]) {
              n1->conn0[j]->t[0] = add_tangent (thisSP->parent->tangents,
                dim, thisSP->t[1]->x, 0);
        } } }
        for (j=0; j<n1->numconn1; j++) {
          if (n1->conn1[j] != thisSP) {
            if (!n1->conn1[j]->t[1]) {
              n1->conn1[j]->t[1] = add_tangent (thisSP->parent->tangents,
                dim, thisSP->t[1]->x, 0);
        } } }

      } else {

        // if this segment is tiny, then just plug it in straight
        // if its half, then just take the mean
        // if its a lot more, then, um, 
        for (i=0; i<dim; i++)
          thisSP->t[1]->x[i] = thisrad*x[i] + (1.-thisrad)*dlnorm[i];
        (void) vec_normalize (thisSP->t[1]->x, thisSP->t[1]->x, dim);
      }
    }
  }

  // if both nodes are ends, then the segment is straight
  if (do_later0 && do_later1) {
    for (i=0; i<dim; i++) thisSP->t[0]->x[i] = dl[i];
    for (i=0; i<dim; i++) thisSP->t[1]->x[i] = thisSP->t[0]->x[i];
    //fprintf(stderr,"    special case 1\n");

  // if just one is an end, copy the other end's tangent, but mirror it
  } else if (do_later0) {

    // use the mirror image of tangent 1
    dotp = 0.;
    for (i=0; i<dim; i++)
      dotp += dlnorm[i] * thisSP->t[1]->x[i];
    for (i=0; i<dim; i++)
      thisSP->t[0]->x[i] = 2.*dlnorm[i]*dotp - thisSP->t[1]->x[i];
    //fprintf(stderr,"    special case 2\n");

  } else if (do_later1) {

    // use the mirror image of tangent 0
    dotp = 0.;
    for (i=0; i<dim; i++)
      dotp += dlnorm[i] * thisSP->t[0]->x[i];
    for (i=0; i<dim; i++)
      thisSP->t[1]->x[i] = 2.*dlnorm[i]*dotp - thisSP->t[0]->x[i];
    //fprintf(stderr,"    special case 3\n");
  }

  return;
}


/*
 * Recursively set all tangents
 */
int set_seg_tangents (seg_group_ptr thisSG, const int dims) {

  seg_ptr curr;

  // there are no children, check all local entries
  curr = thisSG->first;
  while (curr) {
    find_seg_tangents (curr, dims);
    curr = curr->next;
  }

  return (0);
}


/*
 * Perturb all nodes according to Gaussian, scale by mean segment length
 */
int roughen_nodes (node_group_ptr thisNG, int dim, double scale) {

  int i;
  int cnt = 0;
  double thislen;
  double len = 9.9e+9;
  double pert[dim];
  node_ptr curr;

  // check all children
  if (thisNG->child[0]) {
    cnt += roughen_nodes (thisNG->child[0], dim, scale);
    cnt += roughen_nodes (thisNG->child[1], dim, scale);

  } else {
    // check all local nodes
    curr = thisNG->first;
    while (curr) {
      // check all connected segments, find shortest length
      for (i=0; i<curr->numconn0; i++) {
        thislen = seg_length (curr->conn0[i],dim);
        if (thislen < len) len = thislen;
      }
      for (i=0; i<curr->numconn1; i++) {
        thislen = seg_length (curr->conn1[i],dim);
        if (thislen < len) len = thislen;
      }
      (void) perturb_gaussian (pert, scale*len, dim);
      for (i=0; i<dim; i++) curr->x[i] += pert[i];
      curr = curr->next;
    }
  }

  return (cnt);
}


/*
 * Set the block id on each segment, a block is a set of connected segs
 */
int identify_separate_strands (seg_group_ptr thisSG) {

  int i;
  int cnt = 0;
  int blockcnt = 0;
  seg_ptr curr;

  // reset the counters
  thisSG->numblock = 0;
  curr = thisSG->first;
  while (curr) {
    curr->flag = false;
    curr->block = 0;
    curr = curr->next;
  }

  // check all elements
  curr = thisSG->first;
  while (curr) {

    // has this segment been done yet?
    if (!curr->flag) {
      // reset the counter
      cnt = 0;
      // recursively look through the segments, setting the block and flag
      cnt = identify_connected_segments (curr, blockcnt);
      // what happened?
      //fprintf (stderr,"block %d has %d segments\n",blockcnt,cnt);
      // update the ID count
      blockcnt++;
    }
    curr = curr->next;
  }

  thisSG->numblock = blockcnt;

  return (blockcnt);
}


/*
 * Recursively look at all connected segments! Set block
 */
int identify_connected_segments (seg_ptr thisSP, unsigned int blockid) {

  int i;
  int cnt = 1;

  // if this segment has been done, return immediately
  if (thisSP->flag) return (0);

  // set ID
  thisSP->block = blockid;

  // set flag
  thisSP->flag = true;

  // recursively call all connected elems
  for (i=0; i<thisSP->n[0]->numconn0; i++)
    cnt += identify_connected_segments (thisSP->n[0]->conn0[i], blockid);
  for (i=0; i<thisSP->n[0]->numconn1; i++)
    cnt += identify_connected_segments (thisSP->n[0]->conn1[i], blockid);
  for (i=0; i<thisSP->n[1]->numconn0; i++)
    cnt += identify_connected_segments (thisSP->n[1]->conn0[i], blockid);
  for (i=0; i<thisSP->n[1]->numconn1; i++)
    cnt += identify_connected_segments (thisSP->n[1]->conn1[i], blockid);

  return (cnt);
}


/*
 * Given a strand ID, find all nodes in that strand that have only 1 neighbor
 * return one of them
 */
node_ptr find_end_node_of_strand (seg_group_ptr thisSG, const unsigned int blockid) {

  int cnt = 0;
  node_ptr anend = NULL;
  fprintf (stderr,"  looking for end nodes of strand %d...",blockid);

  // check all elements
  node_ptr currn = thisSG->nodes->first;
  while (currn) {
    // identify which group it is in
    unsigned int myblock = -1;
    if (currn->numconn0 > 0) myblock = currn->conn0[0]->block;
    if (currn->numconn1 > 0) myblock = currn->conn1[0]->block;

    // is this node in the desired strand?
    if (myblock == blockid) {
      const int nconnsegs = currn->numconn0 + currn->numconn1;
      if (nconnsegs < 2) {
        anend = currn;
        ++cnt;
      }
    }
    currn = currn->next;
  }
  fprintf (stderr,"found %d\n",cnt);

  return anend;
}


/*
 * Set the flag for each node that is a root of its strand!
 *
 * Make sure each strand has *at least* one root. It may have more.
 *
 * A strand will have one root for each segment that *crosses*
 * the root plane. If no segments cross the root plane, it will pick
 * the node that is closest to the root plane.
 *
 * NO! To make the solution possible, we must only have one root per strand!
 */
int find_root_of_each_strand (seg_group_ptr thisSG, ACTION *args) {

  int i;
  int cnt = 0;
  int axis = args->iarg[0];
  double intercept = args->darg[1];
  double dist;
  unsigned char *found_root;
  double *close_dist;
  node_ptr *close_node;
  seg_ptr curr;

  // set appropriate flags
  (void) set_node_flags (thisSG->nodes, false);

  // create array for each block
  found_root = (unsigned char*) malloc ((int)(thisSG->numblock+1)*sizeof(char));
  close_dist = (double*) malloc ((int)(thisSG->numblock+1) * sizeof(double));
  close_node = (node_ptr*) malloc ((int)(thisSG->numblock+1) * sizeof(node_ptr));

  // initialize those arrays
  for (i=0; i<thisSG->numblock+1; i++) {
    found_root[i] = false;
    close_dist[i] = 9.9e+9;
    close_node[i] = NULL;
  }

  // check all elements
  curr = thisSG->first;
  while (curr) {

    // have we found a root for this strand/block?
    if (!found_root[curr->block]) {

      // does this segment intersect the root plane?
      if ((curr->n[0]->x[axis] - intercept)*
          (curr->n[1]->x[axis] - intercept) < 0.) {

        found_root[curr->block] = true;

        // one of the two nodes will be the root
        if (fabs(curr->n[0]->x[axis] - intercept) >
            fabs(curr->n[1]->x[axis] - intercept)) {
          curr->n[0]->flag = true;
        } else {
          curr->n[1]->flag = true;
        }
        cnt++;

      // if not, how close to the root plane does it get?
      } else {

        dist = fabs(curr->n[0]->x[axis] - intercept);
        if (dist < close_dist[curr->block]) {
          close_dist[curr->block] = dist;
          close_node[curr->block] = curr->n[0];
        }

        dist = fabs(curr->n[1]->x[axis] - intercept);
        if (dist < close_dist[curr->block]) {
          close_dist[curr->block] = dist;
          close_node[curr->block] = curr->n[1];
        }

      }
    }
    curr = curr->next;
  }

  // check for correctness of results
  for (i=1; i<thisSG->numblock+1; i++) {
    // if there's no intersection, there's no root yet
    if (!found_root[i]) {
      if (close_node[i]) {
        close_node[i]->flag = true;
        //fprintf (stderr,"seg %d closest node at %g\n",i,close_dist[i]);
        cnt++;
      } else {
        fprintf (stderr,"ERROR (f_r_o_e_s): no closest node found\n");
        exit(0);
      }
    }
  }

  return (cnt);
}


/*
 * Remove any tip nodes and their connected segments
 */
int prune_once (seg_group_ptr thisSG) {

  int cnt = 0;
  int n0conn,n1conn;
  seg_ptr curr;

  // make sure we don't just zipper delete a long strand
  (void) set_node_flags (thisSG->nodes, false);

  // search structures for nodes with only one connected segment
  curr = thisSG->first;
  while (curr) {

    // if either node has been part of a merge/delete, skip this segment
    if (!curr->n[0]->flag && !curr->n[1]->flag) {

      // how many connected segments does it have
      n0conn = curr->n[0]->numconn0 + curr->n[0]->numconn1;
      n1conn = curr->n[1]->numconn0 + curr->n[1]->numconn1;

      // if one node has only one segment, knock it out
      if (n0conn == 1) {
        // remove that segment and then remove that node
        delete_segment (curr->n[1], curr->n[0], curr, thisSG->dim);
        cnt++;
      } else if (n1conn == 1) {
        // remove that segment and then remove that node
        delete_segment (curr->n[0], curr->n[1], curr, thisSG->dim);
        cnt++;
      }
    }

    curr = curr->next;
  }

  // now prune all unattached nodes

  return (cnt);
}


/*
 * Set as many segment radii as possible
 *
 * Only works for dimensions 2 and 3 (or first 3 dimensions of d>3)
 *
 * A node flag of true means that that node is a root!
 */
int set_treelike_radii_3d (seg_group_ptr thisSG, ACTION *args) {

  int i,j,n0conn,n1conn,rindex,tindex;
  int cnt = 0;
  double lensq,armsq,rad,thisforce,dl[3];
  double **moment;
  double **force;
  double stress = args->darg[0];
  int axis = args->iarg[0];
  seg_ptr curr;
  node_ptr rootnode,tipnode;

  // allocate space for a moment on each node
  int np1 = thisSG->nodes->num+1;
  moment = (double**) malloc ((int)np1 * sizeof(double*));
  for (i=0; i<np1; i++)
    moment[i] = (double*) malloc ((int)thisSG->dim * sizeof(double));
  force = (double**) malloc ((int)np1 * sizeof(double*));
  for (i=0; i<np1; i++)
    force[i] = (double*) malloc ((int)thisSG->dim * sizeof(double));

  // zero those moments (go to n+1 because of 1-indexing)
  for (i=0; i<np1; i++) {
    for (j=0; j<thisSG->dim; j++) {
      moment[i][j] = 0.;
      force[i][j] = 0.;
    }
  }

  // check all elements
  curr = thisSG->first;
  while (curr) {

    // has this segment been done yet?
    if (!curr->flag) {

      // how many connected segments does it have
      n0conn = curr->n[0]->numconn0 + curr->n[0]->numconn1;
      n1conn = curr->n[1]->numconn0 + curr->n[1]->numconn1;
      //fprintf (stderr,"seg %d numconn %d\n",curr->index,numconn);

      // it's independent
      if (n0conn == 1 && n1conn == 1) {

        // it doesn't matter which end is the root
        // segment length
        lensq = seg_length_sqrd (curr, thisSG->dim);
        // moment arm
        armsq = lensq - pow(curr->n[1]->x[axis] - curr->n[0]->x[axis],2);
        // here's the calculation for a constant-radius segment
        //   with no tip force or moment
        rad = 2.*sqrt(lensq*armsq)/stress;
        curr->r[0] = add_radius (thisSG->radii, rad, 0);
        curr->r[1] = add_radius (thisSG->radii, rad, 0);
        //fprintf (stderr,"  independent, has radius %g\n",rad);

      // it's an end, do it
      } else if (n0conn == 1 || n1conn == 1) {

        // define our directions
        if (n1conn == 1) {
          rootnode = curr->n[0];
          rindex = curr->n[0]->index;
          tipnode = curr->n[1];
          tindex = curr->n[1]->index;
        } else {
          rootnode = curr->n[1];
          rindex = curr->n[1]->index;
          tipnode = curr->n[0];
          tindex = curr->n[0]->index;
        }

        // segment length
        for (j=0; j<3; j++) dl[j] = tipnode->x[j] - rootnode->x[j];
        lensq = dl[0]*dl[0] + dl[1]*dl[1] + dl[2]*dl[2];
        // moment arm
        armsq = lensq - pow(curr->n[1]->x[axis] - curr->n[0]->x[axis],2);
        // here's the calculation for a constant-radius segment
        //   with no tip force or moment
        rad = 2.*sqrt(lensq*armsq)/stress;
        curr->r[0] = add_radius (thisSG->radii, rad, 0);
        curr->r[1] = add_radius (thisSG->radii, rad, 0);
        fprintf (stderr,"len %g, arm %g\n",sqrt(lensq),sqrt(armsq));
        fprintf (stderr,"  tip has radius %g\n",rad);

        // set the force and moment for this tip segment
        force[rindex][axis] = M_PI*sqrt(lensq)*rad*rad;
        vec_cross (dl, force[rindex], moment[rindex]);
        fprintf (stderr,"  force %g\n",force[rindex][2]);
        fprintf (stderr,"  moment %g %g\n",moment[rindex][0],moment[rindex][1]);

        // carry this force all the way down to root?

      // it's not an end, check all connected segs for doneness
      } else {

        // UNFINISHED!!!
        // only continue if ONE connected adjacent seg is not done

        // define our directions
        if (n1conn == 1) {
          rootnode = curr->n[0];
          rindex = curr->n[0]->index;
          tipnode = curr->n[1];
          tindex = curr->n[1]->index;
        } else {
          rootnode = curr->n[1];
          rindex = curr->n[1]->index;
          tipnode = curr->n[0];
          tindex = curr->n[0]->index;
        }

        // carry the force from the tipward segments
        //for (j=0; j<3; j++)
          //force[rootnode->index][j] += force[tipnode->index][j];
      }

      // for testing, force all segments to be done
      curr->flag = true;
    }
    curr = curr->next;
  }

  return (cnt);
}


/*
 * Translate all nodes by a fixed vector
 */
int translate_nodes (node_group_ptr thisNG, int dim, double *vec) {

  int i;
  int cnt = 0;
  node_ptr curr;

  // check all children
  if (thisNG->child[0]) {
    cnt += translate_nodes (thisNG->child[0], dim, vec);
    cnt += translate_nodes (thisNG->child[1], dim, vec);

  } else {
    // check all local nodes
    curr = thisNG->first;
    while (curr) {
      for (i=0; i<dim; i++) curr->x[i] += vec[i];
      curr = curr->next;
    }
  }

  return (cnt);
}


/*
 * Scale all nodes by a fixed vector
 */
int scale_nodes (node_group_ptr thisNG, int dim, double *vec) {

  int i;
  int cnt = 0;
  node_ptr curr;

  // check all children
  if (thisNG->child[0]) {
    cnt += scale_nodes (thisNG->child[0], dim, vec);
    cnt += scale_nodes (thisNG->child[1], dim, vec);

  } else {
    // check all local nodes
    curr = thisNG->first;
    while (curr) {
      for (i=0; i<dim; i++) curr->x[i] *= vec[i];
      curr = curr->next;
    }
  }

  // finally, scale the min/max bounds
  for (i=0; i<dim; i++) thisNG->min[i] *= vec[i];
  for (i=0; i<dim; i++) thisNG->max[i] *= vec[i];

  return (cnt);
}


/*
 * Scale all radii by either a fixed number, or interpolate along a range of two numbers
 */
int scale_radii (rad_group_ptr thisRG, double *vec) {

  int i;
  int cnt = 0;
  rad_ptr curr;
  static int setScaling = false;
  static double c1, c2;

  // have we determined the scaling yet?
  if (!setScaling) {
    // how many numbers were given for scaling?
    if (vec[0] == DBLFLAG) {
      // no numbers given, do not scale radii
      c1 = 0.0;
      c2 = 1.0;
    } else if (vec[1] == DBLFLAG) {
      // one number given, uniformly scale all radii
      c1 = 0.0;
      c2 = vec[0];
    } else {
      // two numbers given, scale radii linearly along range
      c2 = (vec[1] - vec[0]) / (thisRG->max - thisRG->min);
      c1 = vec[0] - c2*thisRG->min;
    }
    //thisRG->min *= vec[0];
    //thisRG->max *= vec[0];
    setScaling = true;
  }

  // check all children
  if (thisRG->child[0]) {
    cnt += scale_radii (thisRG->child[0], vec);
    cnt += scale_radii (thisRG->child[1], vec);

  } else {
    // check all local nodes
    curr = thisRG->first;
    while (curr) {
      //curr->r *= vec[0];
      curr->r = c1 + curr->r*c2;
      curr = curr->next;
    }
  }

  // finally, scale the min/max bounds
  thisRG->min = c1 + thisRG->min*c2;
  thisRG->max = c1 + thisRG->max*c2;

  return (cnt);
}


/*
 * Scale all nodes by a fixed vector
 */
int split_and_write_rad (seg_group_ptr thisSG, char root[MAXSTR]) {

  int i;
  int long_axis = -1;
  double long_length = -1.;
  double this_length = -1.;
  char leftroot[MAXSTR],rightroot[MAXSTR];
  FILE *outptr = stdout;

  // two child segment groups
  seg_group_ptr left,right;

  //fprintf(stderr,"splitting (%d)\n",thisSG->num);
  fprintf(stderr,"splitting (%s)\n",root);

  // is this group small enough to just print?
  // construct the filename
  //if (true) {
  if (thisSG->num < 100000) {
    sprintf(leftroot,"%s.rad",root);
    outptr = fopen(leftroot,"w");
    write_rad(outptr,thisSG);
    fclose(outptr);
    // don't delete the structure here because it may be the original segs
    return(-1);
  }

  // if we're still here, then we need to split this seg group into two child seg groups

  // prepare both groups
  left = make_seg_structure (thisSG->dim);
  right = make_seg_structure (thisSG->dim);

  //dump_stats(thisSG);

  // find the longest edge
  long_length = -1.;
  for (i=0; i<thisSG->dim; i++) {
    this_length = pow(thisSG->nodes->max[i]-thisSG->nodes->min[i],2);
    if (this_length > long_length) {
      long_length = this_length;
      long_axis = i;
    }
  }
  long_length = thisSG->nodes->min[long_axis] + 0.5*sqrt(long_length);

  fprintf(stderr,"  axis %d, split position %g\n",long_axis,long_length);

  // recurse over all segments in "this" and add them to left, right, or both
  split_into_two_seg_groups (thisSG->nodes, long_axis, long_length, left, right);
  fprintf(stderr,"    left segments (%d)\n",left->num);
  fprintf(stderr,"    right segments (%d)\n",right->num);
  //exit(0);

  // now, call this routine on both
  // and remove the structure when we're done with it!
  sprintf(leftroot,"%sl",root);
  (void) split_and_write_rad(left, leftroot);
  (void) free_segs (left);

  sprintf(rightroot,"%sr",root);
  (void) split_and_write_rad(right, rightroot);
  (void) free_segs (right);

  return(0);
}


/*
 * Put all segments, nodes, and tangents into either left, right, or both;
 * cnt is number of segments put on the left.
 */
int split_into_two_seg_groups (node_group_ptr thisNG, int split_dim, double split_val,
                               seg_group_ptr left, seg_group_ptr right) {

  int i;
  int cnt = 0;
  node_ptr curr,n0,n1;
  seg_ptr this_seg,newseg;
  rad_ptr r0,r1;
  //tan_ptr t0,t1;
  enum { leftonly,
         rightonly,
         split } which_side;


  // check all children
  if (thisNG->child[0]) {
    cnt += split_into_two_seg_groups (thisNG->child[0], split_dim, split_val, left, right);
    cnt += split_into_two_seg_groups (thisNG->child[1], split_dim, split_val, left, right);

  } else {
    // check all local nodes
    curr = thisNG->first;
    while (curr) {

      // loop over all connected segments for which this is node 0
      for (i=0; i<curr->numconn0; i++) {
        this_seg = curr->conn0[i];

        // figure out which side to put it on
        // for now, make it totally simple
        if (curr->x[split_dim] < split_val) {
          which_side = leftonly;
        } else {
          which_side = rightonly;
        }

        // now, move it there
        if (which_side == leftonly) {
          n0 = add_node (left->nodes, left->dim, curr->x, 0);
          n1 = add_node (left->nodes, left->dim, this_seg->n[1]->x, 0);
          r0 = add_radius (left->radii, this_seg->r[0]->r, 0);
          r1 = add_radius (left->radii, this_seg->r[1]->r, 0);
          //t0 = add_tangent (left->tangents, this_seg->t[0]->x, 0);
          //t1 = add_tangent (left->tangents, this_seg->t[1]->x, 0);
          newseg = add_segment (left, n0, n1);
          if (r0) newseg->r[0] = r0;
          if (r1) newseg->r[1] = r1;
          //if (t0) newseg->t[0] = t0;
          //if (t1) newseg->t[1] = t1;
        } else if (which_side == rightonly) {
          n0 = add_node (right->nodes, right->dim, curr->x, 0);
          n1 = add_node (right->nodes, right->dim, this_seg->n[1]->x, 0);
          r0 = add_radius (right->radii, this_seg->r[0]->r, 0);
          r1 = add_radius (right->radii, this_seg->r[1]->r, 0);
          newseg = add_segment (right, n0, n1);
          if (r0) newseg->r[0] = r0;
          if (r1) newseg->r[1] = r1;
        } else {
          fprintf(stderr,"ERROR (split_into_two_seg_groups): cannot split\n");
          exit(0);
        }

      }

      curr = curr->next;
    }
  }

  return (cnt);
}



// ==========================================================================
// input/output
//

/*
 * Read a .seg file (syntax very much like .obj)
 */
int read_seg (char *infile, seg_group_ptr thisSG, int zero_indexed) {

  int i,j,k,nnode,nrad,ntan,nlines;
  int n0index,n1index,n0rad,n1rad,n0tan,n1tan;
  char onechar,anotherchar;
  char twochar[2],sbuf[128],newval[32],sub[10];
  double tempd[MAXDIM];
  node_ptr newnode,newnode0,newnode1;
  node_ptr *thenodes;
  rad_ptr newrad;
  rad_ptr *therads;
  tan_ptr newtan;
  tan_ptr *thetans;
  seg_ptr newseg;
  FILE *infp;
  long int fpos;

  // before proceeding, prepare the seg_group struct
  thisSG->dim = 3;		// provide a default
  thisSG->radius = -1.;		// negative == "unset"
  thisSG->num = 0;		// gotta start somewhere

  // first, scan the file and record the dimensionality and bounds =======

  // open the file for reading
  infp = fopen (infile,"r");
  if (infp==NULL) {
    fprintf (stderr,"  Could not open input file %s\n",infile);
    fprintf (stderr,"  Exiting\n");
    exit (0);
  }
  fprintf (stderr,"  Opening file %s\n",infile);
  fprintf (stderr,"  Prescanning...");
  fflush (stderr);

  // read the first character after the newline
  nnode = 0;
  nrad = 0;
  ntan = 0;
  nlines = 0;
  while (fread (&onechar,sizeof(char),1,infp) == 1) {

    if (onechar == 'v') {
      // this is some form of vertex information
      fread (&anotherchar,sizeof(char),1,infp);

      if (isspace (anotherchar)) {
        // this is an actual vertex
        // count it
        nnode++;
      } else if (anotherchar == 'r') {
        // this is a radius
        // count it
        nrad++;
      } else if (anotherchar == 'n') {
        // this is a tangent
        // count it
        ntan++;
      } else {
        // skip it
      }

    } else if (onechar == 'd') {
      // the dimensionality
      fscanf (infp,"%s",newval);
      thisSG->dim = (char)atoi(newval);
      if (thisSG->dim > MAXDIM) {
        fprintf (stderr,"ERROR (read_seg): input file's dim is %d, but\n",
                 thisSG->dim);
        fprintf (stderr,"  the maximum spatial dimensions allowed is %d.\n",
                 (MAXDIM));
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }
      //fprintf (stderr,"set dimensionality to %d\n",(int)thisSG->dim);

    } else if (onechar == 'g') {
      // this is some form of global variable
      fread (&anotherchar,sizeof(char),1,infp);

      if (anotherchar == 'r') {
        // the global radius
        fscanf (infp,"%s",newval);
        thisSG->radius = atof(newval);
        //fprintf (stderr,"set radius to %g\n",thisSG->radius);
      } else {
        // nothing important
      }

    } else {
      // if its not identifiable, skip it
    }

    // finish reading the line and the CR
    fscanf (infp,"%[^\n]",sbuf);	// read comment beyond '#'
    fscanf (infp,"%[\n]",twochar);	// read newline

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // close it
  fclose(infp);
  fprintf(stderr,"%d nodes",nnode);
  if (nrad > 0) fprintf(stderr,", %d radii",nrad);
  if (ntan > 0) fprintf(stderr,", %d tangents",ntan);
  fprintf(stderr,"\n");
  fflush(stderr);

  // malloc space for an array of pointers
  thenodes = (node_ptr*) malloc ((nnode+1) * sizeof(node_ptr));
  therads = (rad_ptr*) malloc ((nrad+1) * sizeof(rad_ptr));
  thetans = (tan_ptr*) malloc ((ntan+1) * sizeof(tan_ptr));

  // then, reopen the file and read in the vertexes and segments =========

  // open the file for reading
  infp = fopen (infile,"r");
  if (infp==NULL) {
    fprintf (stderr,"  Could not open input file %s\n",infile);
    fprintf (stderr,"  Exiting\n");
    exit (0);
  }
  //fprintf (stderr,"  Opening file %s\n",infile);
  fprintf (stderr,"  Reading...");
  fflush (stderr);

  // read the first character after the newline
  nnode = 0;
  nrad = 0;
  nlines = 0;
  while (fread (&onechar,sizeof(char),1,infp) == 1) {

    // split on first character
    if (onechar == '#') {
      // read a comment line

    } else if (onechar == 'v') {
      // this is some form of vertex information
      fread (&anotherchar,sizeof(char),1,infp);

      if (isspace(anotherchar)) {
        // this is an actual vertex
        //fprintf (stderr,"adding node\n"); fflush (stderr);
        for (i=0; i<thisSG->dim; i++) {
          fpos = ftell (infp);
          fscanf (infp,"%s",newval);
          // check this value for \n or non-math character!
          //fprintf (stderr,"    (%s)\n",newval); fflush (stderr);
          if (newval[0] < 43 || newval[0] == 47 | newval[0] > 57) {
            // if it is, fill the rest of the data with zeros
            for (j=i; j<thisSG->dim; j++) tempd[j] = 0.;
            // rewind the file
            fseek (infp, fpos, SEEK_SET);
            // and break out of this loop
            break;
          } else {
            tempd[i] = atof(newval);
          }
        }
        newnode = add_node (thisSG->nodes, thisSG->dim, tempd, 0);
        //for (i=0; i<thisSG->dim; i++) newnode->x[i] = tempd[i];
        //fprintf (stderr," node index %ld\n",newnode->index); fflush (stderr);
        //fprintf (stderr," node loc %g %g\n",newnode->x[0],newnode->x[1]); fflush (stderr);
        // put it in the list
        //thenodes[thisSG->nodes->num] = newnode;
        thenodes[++nnode] = newnode;

      } else if (anotherchar == 'r') {
        // this is a radius
        //fprintf (stderr,"adding radius\n"); fflush (stderr);
        fscanf (infp,"%s",newval);
        newrad = add_radius (thisSG->radii, atof(newval), 0);
        therads[++nrad] = newrad;

      } else if (anotherchar == 'n') {
        // this is a tangent
        for (i=0; i<thisSG->dim; i++) {
          fscanf (infp,"%s",newval);
          tempd[i] = atof(newval);
        }
        //newtan = add_tangent (thisSG->tangents, atof(newval), 0);
        //thetans[++ntan] = newtan;

      } else {
        // skip it
      }

    } else if (onechar == 's') {
      // this is a segment
      //fprintf (stderr,"adding segment\n"); fflush (stderr);
      n0index = -1;
      n1index = -1;
      n0rad = -1;
      n1rad = -1;
      n0tan = -1;
      n1tan = -1;

      // read the first set of indices
      fscanf (infp,"%s",newval);

      // look for a slash
      for (i=0; i<strlen(newval); i++)
        if (newval[i] == '/' || newval[i] == '\0') break;
      strncpy(sub,newval,i);
      sub[i] = '\0';
      if (i>0) n0index = atoi(sub);
      //fprintf (stderr,"i %d, sub (%s) (%d)\n",i,sub,atoi(sub));
      // should we continue?
      if (newval[i] == '/') {
        for (j=i+1; j<strlen(newval); j++)
          if (newval[j] == '/' || newval[j] == '\0') break;
        strncpy(sub,newval+i+1,j);
        sub[j-i-1] = '\0';
        if (j-i > 1) n0rad = atoi(sub);
        //fprintf (stderr,"i %d, j %d, sub (%s) (%d)\n",i,j,sub,atoi(sub));
        // should we continue?
        if (newval[j] == '/') {
          for (k=j+1; k<strlen(newval); k++)
            if (newval[k] == '/' || newval[k] == '\0') break;
          strncpy(sub,newval+j+1,k);
          sub[k-j-1] = '\0';
          if (k-j > 1) n0tan = atoi(sub);
          //fprintf (stderr,"j %d, k %d, sub (%s) (%d)\n",j,k,sub,atoi(sub));
        }
      }

      // read the second set of indices
      fscanf (infp,"%s",newval);

      // look for a slash
      for (i=0; i<strlen(newval); i++)
        if (newval[i] == '/' || newval[i] == '\0') break;
      strncpy(sub,newval,i);
      sub[i] = '\0';
      n1index = atoi(sub);
      //fprintf (stderr,"i %d, sub (%s) (%d)\n",i,sub,atoi(sub));
      // should we continue?
      if (newval[i] == '/') {
        for (j=i+1; j<strlen(newval); j++)
          if (newval[j] == '/' || newval[j] == '\0') break;
        strncpy(sub,newval+i+1,j);
        sub[j-i-1] = '\0';
        n1rad = atoi(sub);
        //fprintf (stderr,"i %d, j %d, sub (%s) (%d)\n",i,j,sub,atoi(sub));
        // should we continue?
        if (newval[j] == '/') {
          for (k=j+1; k<strlen(newval); k++)
            if (newval[k] == '/' || newval[k] == '\0') break;
          strncpy(sub,newval+j+1,k);
          sub[k-j-1] = '\0';
          n1tan = atoi(sub);
          //fprintf (stderr,"j %d, k %d, sub (%s) (%d)\n",j,k,sub,atoi(sub));
        }
      }

      //fprintf (stderr,"  indices %d %d %d %d %d %d\n",n0index,n0rad,n0tan,n1index,n1rad,n1tan); fflush (stderr);

      // check for errors
      if (zero_indexed) {
        if (n0index > -1) n0index++;
        if (n1index > -1) n1index++;
        if (n0rad > -1) n0rad++;
        if (n1rad > -1) n1rad++;
        if (n0tan > -1) n0tan++;
        if (n1tan > -1) n1tan++;
      }

      if (n0index < 1 || n1index < 1) {
        fprintf (stderr,"ERROR (read_seg): segment->node indicies incorrect\n");
        fprintf (stderr,"  n0index=%d, n1index=%d\n",n0index,n1index);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }
      if (n0index > nnode || n1index > nnode) {
        fprintf (stderr,"ERROR (read_seg): segment->node indicies incorrect\n");
        fprintf (stderr,"  n0index=%d, n1index=%d\n",n0index,n1index);
        fprintf (stderr,"  only %d nodes have been defined so far\n",nnode);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }
      if (n0rad == 0 || n1rad == 0) {
        fprintf (stderr,"ERROR (read_seg): segment->rad indicies incorrect\n");
        fprintf (stderr,"  n0rad=%d, n1rad=%d\n",n0rad,n1rad);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }
      if (n0rad > nrad || n1rad > nrad) {
        fprintf (stderr,"ERROR (read_seg): segment->rad indicies incorrect\n");
        fprintf (stderr,"  n0rad=%d, n1rad=%d\n",n0rad,n1rad);
        fprintf (stderr,"  only %d rads have been defined so far\n",nrad);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }

      // create the segment
      newseg = add_segment (thisSG, thenodes[n0index], thenodes[n1index]);

      // assign the radius and texture/tangent
      if (n0rad > -1) newseg->r[0] = therads[n0rad];
      if (n1rad > -1) newseg->r[1] = therads[n1rad];
      if (n0tan > -1) newseg->t[0] = thetans[n0tan];
      if (n1tan > -1) newseg->t[1] = thetans[n1tan];

    } else {
      // if its not identifiable, skip it
    }

    // and finish reading
    fscanf (infp,"%[^\n]",sbuf);	// read line up to newline
    fscanf (infp,"%[\n]",twochar);	// read newline

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // close it
  fclose (infp);
  fprintf (stderr,"%d segments\n",thisSG->num);
  fflush (stderr);

  // remove the temporary array
  free (thenodes);
  free (therads);
  free (thetans);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Read a .obj file (and keep the edges)
 */
int read_obj (char *infile, seg_group_ptr thisSG, int zero_indexed) {

  int i,j,k,nnode,nlines;
  int n0index,n1index,n2index;
  char onechar,anotherchar;
  char twochar[2],sbuf[128],newval[32],sub[10];
  double tempd[MAXDIM];
  node_ptr newnode,newnode0,newnode1;
  node_ptr *thenodes;
  seg_ptr newseg;
  FILE *infp;
  long int fpos;

  // before proceeding, prepare the seg_group struct
  thisSG->dim = 3;		// provide a default
  thisSG->radius = -1.;		// negative == "unset"
  thisSG->num = 0;		// gotta start somewhere

  // first, scan the file and record the dimensionality and bounds =======

  // open the file for reading
  infp = fopen (infile,"r");
  if (infp==NULL) {
    fprintf (stderr,"  Could not open input file %s\n",infile);
    fprintf (stderr,"  Exiting\n");
    exit (0);
  }
  fprintf (stderr,"  Opening file %s\n",infile);
  fprintf (stderr,"  Prescanning...");
  fflush (stderr);

  // read the first character after the newline
  nnode = 0;
  nlines = 0;
  while (fread (&onechar,sizeof(char),1,infp) == 1) {

    if (onechar == 'v') {
      // this is some form of vertex information
      fread (&anotherchar,sizeof(char),1,infp);

      if (isspace (anotherchar)) {
        // this is an actual vertex
        // count it
        nnode++;
      } else {
        // skip it
      }

    } else {
      // if its not identifiable, skip it
    }

    // finish reading the line and the CR
    fscanf (infp,"%[^\n]",sbuf);	// read comment beyond '#'
    fscanf (infp,"%[\n]",twochar);	// read newline

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // close it
  fclose(infp);
  fprintf(stderr,"%d nodes",nnode);
  fprintf(stderr,"\n");
  fflush(stderr);

  // malloc space for an array of pointers
  thenodes = (node_ptr*) malloc ((nnode+1) * sizeof(node_ptr));

  // then, reopen the file and read in the vertexes and segments =========

  // open the file for reading
  infp = fopen (infile,"r");
  if (infp==NULL) {
    fprintf (stderr,"  Could not open input file %s\n",infile);
    fprintf (stderr,"  Exiting\n");
    exit (0);
  }
  //fprintf (stderr,"  Opening file %s\n",infile);
  fprintf (stderr,"  Reading...");
  fflush (stderr);

  // read the first character after the newline
  nnode = 0;
  nlines = 0;
  while (fread (&onechar,sizeof(char),1,infp) == 1) {

    // split on first character
    if (onechar == '#') {
      // read a comment line

    } else if (onechar == 'v') {
      // this is some form of vertex information
      fread (&anotherchar,sizeof(char),1,infp);

      if (isspace(anotherchar)) {
        // this is an actual vertex
        //fprintf (stderr,"adding node\n"); fflush (stderr);
        for (i=0; i<thisSG->dim; i++) {
          fpos = ftell (infp);
          fscanf (infp,"%s",newval);
          // check this value for \n or non-math character!
          //fprintf (stderr,"    (%s)\n",newval); fflush (stderr);
          if (newval[0] < 43 || newval[0] == 47 | newval[0] > 57) {
            // if it is, fill the rest of the data with zeros
            for (j=i; j<thisSG->dim; j++) tempd[j] = 0.;
            // rewind the file
            fseek (infp, fpos, SEEK_SET);
            // and break out of this loop
            break;
          } else {
            tempd[i] = atof(newval);
          }
        }
        newnode = add_node (thisSG->nodes, thisSG->dim, tempd, 0);
        //for (i=0; i<thisSG->dim; i++) newnode->x[i] = tempd[i];
        //fprintf (stderr," node index %ld\n",newnode->index); fflush (stderr);
        //fprintf (stderr," node loc %g %g\n",newnode->x[0],newnode->x[1]); fflush (stderr);
        // put it in the list
        //thenodes[thisSG->nodes->num] = newnode;
        thenodes[++nnode] = newnode;

      } else {
        // skip it
      }

    } else if (onechar == 'f') {
      // this is a segment
      //fprintf (stderr,"adding segment\n"); fflush (stderr);
      n0index = -1;
      n1index = -1;
      n2index = -1;

      // read the first set of indices
      fscanf (infp,"%s",newval);

      // look for a slash
      for (i=0; i<strlen(newval); i++)
        if (newval[i] == '/' || newval[i] == '\0') break;
      strncpy(sub,newval,i);
      sub[i] = '\0';
      if (i>0) n0index = atoi(sub);
      //fprintf (stderr,"i %d, sub (%s) (%d)\n",i,sub,atoi(sub));

      // read the second set of indices
      fscanf (infp,"%s",newval);

      // look for a slash
      for (i=0; i<strlen(newval); i++)
        if (newval[i] == '/' || newval[i] == '\0') break;
      strncpy(sub,newval,i);
      sub[i] = '\0';
      n1index = atoi(sub);
      //fprintf (stderr,"i %d, sub (%s) (%d)\n",i,sub,atoi(sub));

      // and the third set of indices
      fscanf (infp,"%s",newval);

      // look for a slash
      for (i=0; i<strlen(newval); i++)
        if (newval[i] == '/' || newval[i] == '\0') break;
      strncpy(sub,newval,i);
      sub[i] = '\0';
      n2index = atoi(sub);
      //fprintf (stderr,"i %d, sub (%s) (%d)\n",i,sub,atoi(sub));

      //fprintf (stderr,"  indices %d %d %d %d %d %d\n",n0index,n0rad,n0tan,n1index,n1rad,n1tan); fflush (stderr);

      // check for errors
      if (zero_indexed) {
        if (n0index > -1) n0index++;
        if (n1index > -1) n1index++;
        if (n2index > -1) n2index++;
      }

      if (n0index < 1 || n1index < 1 || n2index < 1) {
        fprintf (stderr,"ERROR (read_obj): segment->node indicies incorrect\n");
        fprintf (stderr,"  indices=%d, %d, %d\n",n0index,n1index,n2index);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }
      if (n0index > nnode || n1index > nnode || n2index > nnode) {
        fprintf (stderr,"ERROR (read_obj): segment->node indicies incorrect\n");
        fprintf (stderr,"  indices=%d, %d, %d\n",n0index,n1index,n2index);
        fprintf (stderr,"  only %d nodes have been defined so far\n",nnode);
        fprintf (stderr,"  Quitting.\n");
        exit(1);
      }

      // create the segment
      newseg = add_segment (thisSG, thenodes[n0index], thenodes[n1index]);
      newseg = add_segment (thisSG, thenodes[n1index], thenodes[n2index]);
      newseg = add_segment (thisSG, thenodes[n2index], thenodes[n0index]);

    } else {
      // if its not identifiable, skip it
    }

    // and finish reading
    fscanf (infp,"%[^\n]",sbuf);	// read line up to newline
    fscanf (infp,"%[\n]",twochar);	// read newline

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // close it
  fclose (infp);
  fprintf (stderr,"%d segments\n",thisSG->num);
  fflush (stderr);

  // remove the temporary array
  free (thenodes);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Read a .rad file (detects Radiance cones and cylinders only)
 *
 * For now, naive read: will maintain copies of nodes
 */
int read_rad (char *infile, seg_group_ptr thisSG) {

  int i,nlines,nfloats;
  char onechar,anotherchar;
  char twochar[2],sbuf[512],token[14][32];
  double tempd[8];
  node_ptr newnode0,newnode1;
  rad_ptr newrad;
  seg_ptr newseg;
  FILE *infp;

  // before proceeding, prepare the seg_group struct
  thisSG->dim = 3;		// Radiance files are always 3D
  thisSG->radius = -1.;		// negative == "unset"
  thisSG->num = 0;		// gotta start somewhere

  // first, scan the file and record the dimensionality and bounds =======
  //   no need!

  // then, reopen the file and read in the vertexes and segments =========

  // open the file for reading
  infp = fopen (infile,"r");
  if (infp==NULL) {
    fprintf (stderr,"Could not open input file %s\n",infile);
    fprintf (stderr,"Exiting\n");
    exit (0);
  }
  //fprintf (stderr,"Opening file %s\n",infile);
  fprintf (stderr,"Reading...");
  fflush (stderr);

  // read a whole line from the input file
  nlines = 0;
  //while (fscanf(infp,"%[^\n]",sbuf) != EOF) {
  // naw, read only the first three tokens
  while (fscanf(infp,"%s %s %s",token[0],token[1],token[2]) != EOF) {

    // clear out the token buffers!
    //for (i=0; i<3; i++) token[i][0] = '\0';

    // grab the first three words (mat key name)
    //sscanf(sbuf,"%s %s %s",token[0],token[1],token[2]);
    //fscanf(infp,"%s %s %s",token[0],token[1],token[2]);
    // grab the first six words (mat key name nc ni nf)
    //sscanf(sbuf,"%s %s %s %s %s %s",token[0],token[1],token[2],
    //  token[3],token[4],token[5]);
    //sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s",token[0],
    //  token[1],token[2],token[3],token[4],token[5],token[6],token[7],
    //  token[8],token[9],token[10],token[11],token[12],token[13]);
    //fprintf (stdout,"   first token is %s\n",token[0]); fflush(stdout);
    //for (i=0; i<3; i++) fprintf (stdout,"%d %s\n",i,token[i]); fflush(stdout);
//exit(0);

    // split on those three tokens only
    if (token[0][0] == '#' || token[0][0] == '!') {
      // read a comment line, or some other descriptor
      fscanf (infp,"%[^\n]",sbuf);	// read line up to newline
      fscanf (infp,"%[\n]",twochar);	// read up to newline
      // fprintf (stdout,"%s\n",sbuf);	// write comment

    } else if (strncmp(token[1],"cylinder",8) == 0) {
      // read in and process a cylinder

      // read in the next 7 tokens
      //sscanf(sbuf,"%s %s %s %s %s %s %s",token[0],token[1],token[2],
      //  token[3],token[4],token[5],token[6]);
      //sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s",token[0],
      //  token[1],token[2],token[3],token[4],token[5],token[6],token[7],
      //  token[8],token[9],token[10],token[11],token[12]);

      // call subroutine to find 7 float values
      nfloats = read_radiance_floats (infp, tempd);
      if (nfloats != 7) {
        fprintf (stderr,"ERROR (read_rad): too many arguments for\n");
        fprintf (stderr,"  cylinder (%d).\n",nfloats);
        fprintf (stderr,"Quitting.\n");
        exit(1);
      }
      //for (i=0; i<7; i++) fprintf (stdout,"%d %g\n",i,tempd[i]); fflush(stdout);

      // create nodes at the two endpoints
      newnode0 = add_node (thisSG->nodes, thisSG->dim, &tempd[0], 0);
      //for (i=0; i<thisSG->dim; i++) newnode0->x[i] = tempd[i];

      newnode1 = add_node (thisSG->nodes, thisSG->dim, &tempd[3], 0);
      //for (i=0; i<thisSG->dim; i++) newnode1->x[i] = tempd[3+i];

      // create the segment
      newseg = add_segment (thisSG, newnode0, newnode1);

      // and create one radius
      newrad = add_radius (thisSG->radii, tempd[6], 0);
      newseg->r[0] = newrad;
      newseg->r[1] = newrad;

      // and finish reading
      fscanf (infp,"%[^\n]",sbuf);	// read line up to newline
      fscanf (infp,"%[\n]",twochar);	// read newline

    } else if (strncmp(token[1],"cone",4) == 0) {
      // read in and process a cone

      // read in the next 8 tokens
      //sscanf(sbuf,"%s %s %s %s %s %s %s %s",token[0],token[1],token[2],
      //  token[3],token[4],token[5],token[6],token[7]);
      //sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s",token[0],
      //  token[1],token[2],token[3],token[4],token[5],token[6],token[7],
      //  token[8],token[9],token[10],token[11],token[12],token[13]);

      // call subroutine to find 8 float values
      nfloats = read_radiance_floats (infp, tempd);
      if (nfloats != 8) {
        fprintf (stderr,"ERROR (read_rad): too many arguments for\n");
        fprintf (stderr,"  cone (%d).\n",nfloats);
        fprintf (stderr,"Quitting.\n");
        exit(1);
      }
      //for (i=0; i<8; i++) fprintf (stdout,"%d %g\n",i,tempd[i]); fflush(stdout);

      // create nodes at the two endpoints
      newnode0 = add_node (thisSG->nodes, thisSG->dim, &tempd[0], 0);
      //for (i=0; i<thisSG->dim; i++) newnode0->x[i] = tempd[i];

      newnode1 = add_node (thisSG->nodes, thisSG->dim, &tempd[3], 0);
      //for (i=0; i<thisSG->dim; i++) newnode1->x[i] = tempd[3+i];

      // create the segment
      newseg = add_segment (thisSG, newnode0, newnode1);

      // and create two radii
      newrad = add_radius (thisSG->radii, tempd[6], 0);
      newseg->r[0] = newrad;

      newrad = add_radius (thisSG->radii, tempd[7], 0);
      newseg->r[1] = newrad;

      // and finish reading
      fscanf (infp,"%[^\n]",sbuf);	// read line up to newline
      fscanf (infp,"%[\n]",twochar);	// read newline

    } else {
      // ignore anything else
      fscanf (infp,"%[^\n]",sbuf);	// read comment beyond '#'
      fscanf (infp,"%[\n]",twochar);	// read newline
    }

    // write the counter
    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }

    // clear out the token buffers!
    for (i=0; i<3; i++) token[i][0] = '\0';
  }

  // close it
  fclose (infp);
  fprintf (stderr,"%d segments\n",thisSG->num);
  fflush (stderr);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Parse a Radiance file, passing over whitespace, to get to the
 * floating-point parameters
 */
int read_radiance_floats (FILE *infp, double *x) {

  char token[32];
  int num = 0;
  int i;

  // read the character args
  (void) fscanf(infp,"%s",token);
  num = atoi(token);
  for (i=0; i<num; i++) fscanf(infp,"%s",token);
  // read the integer args
  (void) fscanf(infp,"%s",token);
  num = atoi(token);
  for (i=0; i<num; i++) fscanf(infp,"%s",token);
  // read the float args
  (void) fscanf(infp,"%s",token);
  num = atoi(token);
  for (i=0; i<num; i++) {
    fscanf(infp,"%s",token);
    x[i] = atof(token);
  }

  return (num);
}


/*
 * Generic procedure for creating a new segment
 */
seg_ptr add_segment (seg_group_ptr thisSG, node_ptr n0, node_ptr n1) {

  int i = 0;
  seg_ptr newseg = NULL;
  seg_ptr *newconnlist;

  // does this segment already exist?
  // check neibs of first node for second node
  int exists = false;
  for (unsigned int i=0; i<n0->numconn0; ++i) {
    if (n0->conn0[i]->n[1] == n1) exists = true;
  }
  for (unsigned int i=0; i<n0->numconn1; ++i) {
    if (n0->conn1[i]->n[0] == n1) exists = true;
  }
  // we don't need to check the other

  if (exists) return NULL;

  // get memory for the node
  newseg = (SEGMENT*) malloc (sizeof(SEGMENT));

  // add it to the head of the list
  newseg->next = thisSG->first;
  newseg->prev = NULL;
  if (newseg->next) newseg->next->prev = newseg;
  thisSG->first = newseg;
  newseg->parent = thisSG;

  //if (n0) fprintf (stderr, "  node 0 exists, index %d\n",n0->index);
  //else fprintf (stderr, "  node 0 does not exist, why?\n");
  //if (n1) fprintf (stderr, "  node 1 exists, index %d\n",n1->index);
  //else fprintf (stderr, "  node 1 does not exist, why?\n");

  // null some pointers
  newseg->n[0] = n0;
  newseg->n[1] = n1;
  newseg->t[0] = NULL;
  newseg->t[1] = NULL;
  newseg->r[0] = NULL;
  newseg->r[1] = NULL;

  // first node is index 1
  newseg->index = ++thisSG->num;
  //fprintf (stderr,"\ncreated segment %ld\n",newseg->index);

  // set the nodes' segment pointers
  //fprintf (stderr,"  n0 had %d conn0\n",n0->numconn0);
  n0->numconn0++;
  // alloc space for new conn list
  newconnlist = (seg_ptr*) malloc (n0->numconn0 * sizeof(seg_ptr));
  // copy data from old conn list
  for (i=0; i<n0->numconn0-1; i++) {
    newconnlist[i] = n0->conn0[i];
    //fprintf (stderr,"    one is %ld\n",n0->conn0[i]->index);
  }
  newconnlist[n0->numconn0-1] = newseg;
  // delete old list
  free (n0->conn0);
  // copy new list
  n0->conn0 = newconnlist;
  //fprintf (stderr,"    now has %d conn0\n",n0->numconn0);
  //fprintf (stderr,"      new one is %ld\n",n0->conn0[n0->numconn0-1]->index);

  // now do it for node 1
  //fprintf (stderr,"  n1 had %d conn1\n",n1->numconn1);
  n1->numconn1++;
  // alloc space for new conn list
  newconnlist = (seg_ptr*) malloc (n1->numconn1 * sizeof(seg_ptr));
  // copy data from old conn list
  for (i=0; i<n1->numconn1-1; i++) {
    newconnlist[i] = n1->conn1[i];
    //fprintf (stderr,"    one is %ld\n",n1->conn1[i]->index);
  }
  newconnlist[n1->numconn1-1] = newseg;
  // delete old list
  free (n1->conn1);
  // copy new list
  n1->conn1 = newconnlist;

  // return the pointer
  return (newseg);
}


/*
 * Generic procedure for creating a new node in a node group
 *
 * Does not do any tree stuff yet
 */
/*
node_ptr add_node (node_group_ptr thisNG, unsigned char dim) {

  int i = 0;
  node_ptr newnode = NULL;

  // if there no nodes, set the bounds
  if (!thisNG->first) {
    thisNG->min = (double*) malloc ((int)dim * sizeof(double));
    thisNG->max = (double*) malloc ((int)dim * sizeof(double));
    for (i=0; i<dim; i++) thisNG->min[i] = 9.9e+9;
    for (i=0; i<dim; i++) thisNG->max[i] = -9.9e+9;
  }

  // get memory for the node
  newnode = (NODE*) malloc (sizeof(NODE));

  // add it to the head of the list
  newnode->next = thisNG->first;
  newnode->prev = NULL;
  if (newnode->next) newnode->next->prev = newnode;
  thisNG->first = newnode;

  // assign some shit
  newnode->numconn0 = 0;
  newnode->conn0 = NULL;
  newnode->numconn1 = 0;
  newnode->conn1 = NULL;
  newnode->x = (double*) malloc ((int)dim * sizeof(double));

  // first node is index 1
  newnode->index = ++thisNG->num;

  // return the pointer
  return (newnode);
}
*/


/*
 * Generic procedure for creating a new node in a node group
 */
node_ptr add_node (node_group_ptr thisNG, unsigned char dim, double *x,
  int totnode) {

  int i = 0;
  int axis;
  double dist;
  node_ptr newnode = NULL;
  node_ptr curr = NULL;
  node_ptr currnext = NULL;

  // set the index of the new entry, if we use it
  if (totnode == 0) totnode = thisNG->num;

  // first, see if we have this entry already
  if (thisNG->child[0]) {

    // there are children, check in the appropriate one
    if (2.*x[thisNG->axis] <
        thisNG->child[0]->max[thisNG->axis] + thisNG->child[1]->min[thisNG->axis]) {
      newnode = add_node (thisNG->child[0], dim, x, totnode);
      // reset the bounds
      for (i=0; i<dim; i++) {
        if (thisNG->child[0]->min[i] < thisNG->min[i])
          thisNG->min[i] = thisNG->child[0]->min[i];
        if (thisNG->child[0]->max[i] > thisNG->max[i])
          thisNG->max[i] = thisNG->child[0]->max[i];
      }
    } else {
      newnode = add_node (thisNG->child[1], dim, x, totnode);
      // reset the bounds
      for (i=0; i<dim; i++) {
        if (thisNG->child[1]->min[i] < thisNG->min[i])
          thisNG->min[i] = thisNG->child[1]->min[i];
        if (thisNG->child[1]->max[i] > thisNG->max[i])
          thisNG->max[i] = thisNG->child[1]->max[i];
      }
    }

    // reset the count
    thisNG->num = thisNG->child[0]->num + thisNG->child[1]->num;

  } else {

    // there are no children, check all local entries
    curr = thisNG->first;
    while (curr) {
      //fprintf (stderr,"  checking vs %g\n",curr->r);
      dist = 0;
      for (i=0; i<dim; i++) dist += pow(x[i]-curr->x[i],2);
      if (fabs(dist) < EPSILONSQRD) {
        // these are the same node!
        newnode = curr;
        //fprintf (stderr,"  found match %g\n",curr->r);
      }
      curr = curr->next;
    }
  }

  // if we couldn't find a match, make a new entry
  if (!newnode) {

    //fprintf (stderr,"  no match, creating %g\n",rad);

    // get memory for the node
    newnode = (NODE*) malloc (sizeof(NODE));

    // assign some shit
    newnode->numconn0 = 0;
    newnode->conn0 = NULL;
    newnode->numconn1 = 0;
    newnode->conn1 = NULL;
    newnode->x = (double*) malloc ((int)dim * sizeof(double));

    // set the value
    for (i=0; i<dim; i++) newnode->x[i] = x[i];

    // if this is the first node in this group, malloc the space
    if (!thisNG->first) {
      thisNG->axis = -1;
      thisNG->min = (double*) malloc ((int)dim * sizeof(double));
      thisNG->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisNG->min[i] = x[i];
      for (i=0; i<dim; i++) thisNG->max[i] = x[i];
    }

    // add it to the head of the list
    newnode->next = thisNG->first;
    thisNG->first = newnode;
    newnode->prev = NULL;
    if (newnode->next) newnode->next->prev = newnode;
    newnode->parent = thisNG;

    // set the index
    newnode->index = totnode + 1;

    // and increment the counter for this group
    thisNG->num++;

    // reset the bounds
    for (i=0; i<dim; i++) {
      if (x[i] < thisNG->min[i]) thisNG->min[i] = x[i];
      if (x[i] > thisNG->max[i]) thisNG->max[i] = x[i];
    }

    // finally, if we surpassed the bucket size, split this group!
    if (thisNG->num > BUCKET) {

      // make two children
      thisNG->child[0] = (NODE_GROUP*) malloc (sizeof(NODE_GROUP));
      thisNG->child[1] = (NODE_GROUP*) malloc (sizeof(NODE_GROUP));

      // choose which axis to split (use dist as the max size here)
      dist = 0.;
      for (i=0; i<dim; i++)
        if (thisNG->max[i] - thisNG->min[i] > dist) {
          dist = thisNG->max[i] - thisNG->min[i];
          //fprintf (stderr," %d %g\n",i,dist);
          axis = i;
        }
      thisNG->axis = axis;

      // dump a sentence
      //fprintf (stderr,"\noriginal: num=%d, bounds %g %g, split axis %d\n",
      //         thisNG->num,thisNG->min[axis],thisNG->max[axis],axis);

      // children point to parent
      thisNG->child[0]->num = 0;
      thisNG->child[0]->first = NULL;
      thisNG->child[0]->parent = thisNG;
      thisNG->child[0]->child[0] = NULL;
      thisNG->child[0]->child[1] = NULL;
      thisNG->child[0]->axis = -1;
      thisNG->child[0]->min = (double*) malloc ((int)dim * sizeof(double));
      thisNG->child[0]->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisNG->child[0]->min[i] = 9.9e+9;
      for (i=0; i<dim; i++) thisNG->child[0]->max[i] = -9.9e+9;
      thisNG->child[1]->num = 0;
      thisNG->child[1]->first = NULL;
      thisNG->child[1]->parent = thisNG;
      thisNG->child[1]->child[0] = NULL;
      thisNG->child[1]->child[1] = NULL;
      thisNG->child[1]->axis = -1;
      thisNG->child[1]->min = (double*) malloc ((int)dim * sizeof(double));
      thisNG->child[1]->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisNG->child[1]->min[i] = 9.9e+9;
      for (i=0; i<dim; i++) thisNG->child[1]->max[i] = -9.9e+9;

      // split the cell along this edge
      dist = 0.5 * (thisNG->max[axis] + thisNG->min[axis]);

      // move all entries to one of the two children
      curr = thisNG->first;
      while (curr) {
        currnext = curr->next;
        if (curr->x[axis] < dist) {
          curr->next = thisNG->child[0]->first;
          thisNG->child[0]->first = curr;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisNG->child[0]->num++;
          curr->parent = thisNG->child[0];
          for (i=0; i<dim; i++)
            if (curr->x[i] < thisNG->child[0]->min[i])
              thisNG->child[0]->min[i] = curr->x[i];
          for (i=0; i<dim; i++)
            if (curr->x[i] > thisNG->child[0]->max[i])
              thisNG->child[0]->max[i] = curr->x[i];
        } else {
          curr->next = thisNG->child[1]->first;
          thisNG->child[1]->first = curr;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisNG->child[1]->num++;
          curr->parent = thisNG->child[1];
          for (i=0; i<dim; i++)
            if (curr->x[i] < thisNG->child[1]->min[i])
              thisNG->child[1]->min[i] = curr->x[i];
          for (i=0; i<dim; i++)
            if (curr->x[i] > thisNG->child[1]->max[i])
              thisNG->child[1]->max[i] = curr->x[i];
        }
        curr = currnext;
      }

      // blank out the entry list
      thisNG->first = NULL;

      // update counter and min/max
      //fprintf (stderr,"  makes children with num=%d, bounds %g %g\n",
      //         thisNG->child[0]->num,thisNG->child[0]->min[axis],
      //         thisNG->child[0]->max[axis]);
      //fprintf (stderr,"                  and num=%d, bounds %g %g\n",
      //         thisNG->child[1]->num,thisNG->child[1]->min[axis],
      //         thisNG->child[1]->max[axis]);
    }
  }

  //if (newnode) {
  //fprintf (stderr,"Added a node %d/%d at %g %g %g\n",totnode,newnode->index,x[0],x[1],x[2]);
  //fprintf (stderr,"Added node %d\n",newnode->index);
  //} else {
  //  fprintf (stderr,"Did not add a node? %d\n",totnode);
  //}

  // return the pointer
  return (newnode);
}


/*
 * Generic procedure for creating a new radius in a radius group
 */
rad_ptr add_radius (rad_group_ptr thisRG, double rad, int totrad) {

  double midpt;
  rad_ptr newrad = NULL;
  rad_ptr curr = NULL;
  rad_ptr currnext = NULL;

  // set the index of the new entry, if we use it
  if (totrad == 0) totrad = thisRG->num;

  // first, see if we have this entry already
  if (thisRG->child[0]) {

    // there are children, check in the appropriate one
    if (2.*rad < thisRG->child[0]->max + thisRG->child[1]->min) {
      newrad = add_radius (thisRG->child[0], rad, totrad);
      // reset the bounds
      if (thisRG->child[0]->min < thisRG->min) thisRG->min = thisRG->child[0]->min;
      if (thisRG->child[0]->max > thisRG->max) thisRG->max = thisRG->child[0]->max;
    } else {
      newrad = add_radius (thisRG->child[1], rad, totrad);
      // reset the bounds
      if (thisRG->child[1]->min < thisRG->min) thisRG->min = thisRG->child[1]->min;
      if (thisRG->child[1]->max > thisRG->max) thisRG->max = thisRG->child[1]->max;
    }

    // reset the count
    thisRG->num = thisRG->child[0]->num + thisRG->child[1]->num;

  } else {

    // there are no children, check all local entries
    curr = thisRG->first;
    while (curr) {
      //fprintf (stderr,"  checking vs %g\n",curr->r);
      if (fabs(rad-curr->r) < EPSILON) {
        // these are the same node!
        newrad = curr;
        //fprintf (stderr,"  found match %g\n",curr->r);
      }
      curr = curr->next;
    }
  }

  // if we couldn't find a match, make a new entry
  if (!newrad) {

    //fprintf (stderr,"  no match, creating %g\n",rad);

    // get memory for the node
    newrad = (RADIUS*) malloc (sizeof(RADIUS));

    // set the value
    newrad->r = rad;

    // add it to the head of the list
    newrad->next = thisRG->first;
    thisRG->first = newrad;
    if (newrad->next) newrad->next->prev = newrad;
    newrad->prev = NULL;
    newrad->parent = thisRG;

    // set the index
    newrad->index = totrad + 1;

    // and increment the counter for this group
    thisRG->num++;

    // reset the bounds
    if (rad < thisRG->min) thisRG->min = rad;
    if (rad > thisRG->max) thisRG->max = rad;

    // finally, if we surpassed the bucket size, split this group!
    if (thisRG->num > BUCKET) {
      //fprintf (stderr,"\noriginal box with num=%d, bounds %g %g\n",
      //         thisRG->num,thisRG->min,thisRG->max);

      // make two children
      thisRG->child[0] = (RADIUS_GROUP*) malloc (sizeof(RADIUS_GROUP));
      thisRG->child[1] = (RADIUS_GROUP*) malloc (sizeof(RADIUS_GROUP));

      // children point to parent
      thisRG->child[0]->num = 0;
      thisRG->child[0]->first = NULL;
      thisRG->child[0]->parent = thisRG;
      thisRG->child[0]->child[0] = NULL;
      thisRG->child[0]->child[1] = NULL;
      thisRG->child[0]->min = thisRG->min;
      thisRG->child[0]->max = thisRG->min;
      thisRG->child[1]->num = 0;
      thisRG->child[1]->first = NULL;
      thisRG->child[1]->parent = thisRG;
      thisRG->child[1]->child[0] = NULL;
      thisRG->child[1]->child[1] = NULL;
      thisRG->child[1]->min = thisRG->max;
      thisRG->child[1]->max = thisRG->max;

      // split the cell along this edge
      midpt = 0.5 * (thisRG->max + thisRG->min);

      // move all entries to one of the two children
      curr = thisRG->first;
      while (curr) {
        currnext = curr->next;
        if (curr->r < midpt) {
          curr->next = thisRG->child[0]->first;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisRG->child[0]->first = curr;
          thisRG->child[0]->num++;
          curr->parent = thisRG->child[0];
          if (curr->r > thisRG->child[0]->max) thisRG->child[0]->max = curr->r;
        } else {
          curr->next = thisRG->child[1]->first;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisRG->child[1]->first = curr;
          thisRG->child[1]->num++;
          curr->parent = thisRG->child[1];
          if (curr->r < thisRG->child[1]->min) thisRG->child[1]->min = curr->r;
        }
        curr = currnext;
      }

      // blank out the entry list
      thisRG->first = NULL;

      // update counter and min/max
      //fprintf (stderr,"  makes children with num=%d, bounds %g %g\n",
      //         thisRG->child[0]->num,thisRG->child[0]->min,thisRG->child[0]->max);
      //fprintf (stderr,"                  and num=%d, bounds %g %g\n",
      //         thisRG->child[1]->num,thisRG->child[1]->min,thisRG->child[1]->max);
    }
  }

  // return the pointer
  return (newrad);
}


/*
 * Generic procedure for creating a new tangent in a tangent group
 */
tan_ptr add_tangent (tan_group_ptr thisTG, unsigned char dim, double *x,
  int tottan) {

  int i = 0;
  int axis;
  double dist;
  tan_ptr newtan = NULL;
  tan_ptr curr = NULL;
  tan_ptr currnext = NULL;

  if (isnan(x[0])) exit(0);
  if (isnan(x[1])) exit(0);

  // set the index of the new entry, if we use it
  if (tottan == 0) tottan = thisTG->num;

  // first, see if we have this entry already
  if (thisTG->child[0]) {

    // there are children, check in the appropriate one
    if (2.*x[thisTG->axis] <
        thisTG->child[0]->max[thisTG->axis] + thisTG->child[1]->min[thisTG->axis]) {
      newtan = add_tangent (thisTG->child[0], dim, x, tottan);
      // reset the bounds
      for (i=0; i<dim; i++) {
        if (thisTG->child[0]->min[i] < thisTG->min[i])
          thisTG->min[i] = thisTG->child[0]->min[i];
        if (thisTG->child[0]->max[i] > thisTG->max[i])
          thisTG->max[i] = thisTG->child[0]->max[i];
      }
    } else {
      newtan = add_tangent (thisTG->child[1], dim, x, tottan);
      // reset the bounds
      for (i=0; i<dim; i++) {
        if (thisTG->child[1]->min[i] < thisTG->min[i])
          thisTG->min[i] = thisTG->child[1]->min[i];
        if (thisTG->child[1]->max[i] > thisTG->max[i])
          thisTG->max[i] = thisTG->child[1]->max[i];
      }
    }

    // reset the count
    thisTG->num = thisTG->child[0]->num + thisTG->child[1]->num;

  } else {

    // there are no children, check all local entries
    curr = thisTG->first;
    while (curr) {
      //fprintf (stderr,"  checking vs %g\n",curr->r);
      dist = 0;
      for (i=0; i<dim; i++) dist += pow(x[i]-curr->x[i],2);
      if (fabs(dist) < EPSILONSQRD) {
        // these are the same tan!
        newtan = curr;
        //fprintf (stderr,"  found match %g\n",curr->r);
      }
      curr = curr->next;
    }
  }

  // if we couldn't find a match, make a new entry
  if (!newtan) {

    fprintf (stderr,"  no match, creating %g %g\n",x[0],x[1]);

    // get memory for the tan
    newtan = (TANGENT*) malloc (sizeof(TANGENT));

    // assign some shit
    newtan->x = (double*) malloc ((int)dim * sizeof(double));

    // set the value
    for (i=0; i<dim; i++) newtan->x[i] = x[i];

    // if thisTG is the first tan in this group, malloc the space
    if (!thisTG->first) {
      thisTG->axis = -1;
      thisTG->min = (double*) malloc ((int)dim * sizeof(double));
      thisTG->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisTG->min[i] = x[i];
      for (i=0; i<dim; i++) thisTG->max[i] = x[i];
    }

    // add it to the head of the list
    newtan->next = thisTG->first;
    thisTG->first = newtan;
    newtan->prev = NULL;
    if (newtan->next) newtan->next->prev = newtan;
    newtan->parent = thisTG;

    // set the index
    newtan->index = tottan + 1;

    // and increment the counter for this group
    thisTG->num++;

    // reset the bounds
    for (i=0; i<dim; i++) {
      if (x[i] < thisTG->min[i]) thisTG->min[i] = x[i];
      if (x[i] > thisTG->max[i]) thisTG->max[i] = x[i];
    }

    // finally, if we surpassed the bucket size, split this group!
    if (thisTG->num > BUCKET) {

      // make two children
      thisTG->child[0] = (TANGENT_GROUP*) malloc (sizeof(TANGENT_GROUP));
      thisTG->child[1] = (TANGENT_GROUP*) malloc (sizeof(TANGENT_GROUP));

      // choose which axis to split (use dist as the max size here)
      dist = 0.;
      for (i=0; i<dim; i++) {
        fprintf (stderr," %d %g %g\n",i,thisTG->min[i],thisTG->max[i]);
        if (thisTG->max[i] - thisTG->min[i] > dist) {
          dist = thisTG->max[i] - thisTG->min[i];
          fprintf (stderr,"   %d %g\n",i,dist);
          axis = i;
        }
      }
      thisTG->axis = axis;

      // dump a sentence
      //fprintf (stderr,"\noriginal: num=%d, bounds %g %g, split axis %d\n",
      //         thisTG->num,thisTG->min[axis],thisTG->max[axis],axis);

      // children point to parent
      thisTG->child[0]->num = 0;
      thisTG->child[0]->first = NULL;
      thisTG->child[0]->parent = thisTG;
      thisTG->child[0]->child[0] = NULL;
      thisTG->child[0]->child[1] = NULL;
      thisTG->child[0]->axis = -1;
      thisTG->child[0]->min = (double*) malloc ((int)dim * sizeof(double));
      thisTG->child[0]->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisTG->child[0]->min[i] = 9.9e+9;
      for (i=0; i<dim; i++) thisTG->child[0]->max[i] = -9.9e+9;
      thisTG->child[1]->num = 0;
      thisTG->child[1]->first = NULL;
      thisTG->child[1]->parent = thisTG;
      thisTG->child[1]->child[0] = NULL;
      thisTG->child[1]->child[1] = NULL;
      thisTG->child[1]->axis = -1;
      thisTG->child[1]->min = (double*) malloc ((int)dim * sizeof(double));
      thisTG->child[1]->max = (double*) malloc ((int)dim * sizeof(double));
      for (i=0; i<dim; i++) thisTG->child[1]->min[i] = 9.9e+9;
      for (i=0; i<dim; i++) thisTG->child[1]->max[i] = -9.9e+9;

      // split the cell along this edge
      dist = 0.5 * (thisTG->max[axis] + thisTG->min[axis]);

      // move all entries to one of the two children
      curr = thisTG->first;
      while (curr) {
        currnext = curr->next;
        if (curr->x[axis] < dist) {
          curr->next = thisTG->child[0]->first;
          thisTG->child[0]->first = curr;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisTG->child[0]->num++;
          curr->parent = thisTG->child[0];
          for (i=0; i<dim; i++)
            if (curr->x[i] < thisTG->child[0]->min[i])
              thisTG->child[0]->min[i] = curr->x[i];
          for (i=0; i<dim; i++)
            if (curr->x[i] > thisTG->child[0]->max[i])
              thisTG->child[0]->max[i] = curr->x[i];
        } else {
          curr->next = thisTG->child[1]->first;
          thisTG->child[1]->first = curr;
          if (curr->next) curr->next->prev = curr;
          curr->prev = NULL;
          thisTG->child[1]->num++;
          curr->parent = thisTG->child[1];
          for (i=0; i<dim; i++)
            if (curr->x[i] < thisTG->child[1]->min[i])
              thisTG->child[1]->min[i] = curr->x[i];
          for (i=0; i<dim; i++)
            if (curr->x[i] > thisTG->child[1]->max[i])
              thisTG->child[1]->max[i] = curr->x[i];
        }
        curr = currnext;
      }

      // blank out the entry list
      thisTG->first = NULL;

      // update counter and min/max
      //fprintf (stderr,"  makes children with num=%d, bounds %g %g\n",
      //         thisTG->child[0]->num,thisTG->child[0]->min[axis],
      //         thisTG->child[0]->max[axis]);
      //fprintf (stderr,"                  and num=%d, bounds %g %g\n",
      //         thisTG->child[1]->num,thisTG->child[1]->min[axis],
      //         thisTG->child[1]->max[axis]);
    }
  }

  //if (newtan) {
  //fprintf (stderr,"Added a tan %d/%d at %g %g %g\n",tottan,newtan->index,x[0],x[1],x[2]);
  //fprintf (stderr,"Added tan %d\n",newtan->index);
  //} else {
  //  fprintf (stderr,"Did not add a tan? %d\n",tottan);
  //}

  // return the pointer
  return (newtan);
}


/*
 * Search for the closest node to a point, or return NULL if there are
 * no nodes within the value of dist given in the initial call.
 */
node_ptr find_closest_node (node_group_ptr thisNG, unsigned char dim, double *x,
  node_ptr closenode, double *dist) {

  int i;
  double thisdistsq;
  double distsq = pow(*dist,2);
  //node_ptr closenode = NULL;
  node_ptr curr = NULL;

  // first, can any nodes in this group be close enough?
  if (thisNG->first != NULL || thisNG->child[0] != NULL) {
    thisdistsq = 0;
    for (i=0; i<dim; i++) thisdistsq +=
      pow( fmax( 0., fmax(x[i]-thisNG->max[i], thisNG->min[i]-x[i]) ) ,2);
    // if no nodes can possibly be closer, return
    if (thisdistsq > distsq) return (closenode);
  }

  // if this group could contain a closer node, check the contents
  if (thisNG->child[0]) {

    // there are children, check both
    closenode = find_closest_node (thisNG->child[0], dim, x, closenode, dist);
    closenode = find_closest_node (thisNG->child[1], dim, x, closenode, dist);

  } else {

    // there are no children, check all local entries
    curr = thisNG->first;
    while (curr) {
      thisdistsq = 0;
      for (i=0; i<dim; i++) thisdistsq += pow(x[i]-curr->x[i],2);
      if (thisdistsq < distsq) {
        // this is the new closest node!
        closenode = curr;
        distsq = thisdistsq;
      }
      curr = curr->next;
    }

    // reset the distance
    *dist = sqrt(distsq);
  }

  // return the pointer
  return (closenode);
}


/*
 * Write a .seg file (syntax very much like .obj)
 */
int write_seg (FILE *out, seg_group_ptr thisSG, int argc, char **argv) {

  int i,j;
  int nnode = 0;
  int nrad = 0;
  int ntan = 0;
  int nlines = 0;
  node_ptr currn;
  rad_ptr currr;
  tan_ptr currt;
  seg_ptr curr;

  fprintf (stderr,"  Writing .seg file");
  fflush (stderr);

  // set all flags to true (true==not been printed yet)
  (void) set_radius_flags (thisSG->radii, true);
  (void) set_node_flags (thisSG->nodes, true);

  // for this first take, blindly dump all segments and nodes

  fprintf (out,"# %s\n",VERSION);
  fprintf (out,"#");
  for (i=0; i<argc; i++) fprintf (out," %s",argv[i]);
  fprintf (out,"\n");

  // dimension and global radius
  fprintf (out,"d %d\n",thisSG->dim);
  fprintf (out,"gr %g\n",thisSG->radius);

  // first, the nodes
/*
  currn = thisSG->nodes->first;
  while (currn) {
    fprintf (out,"v");
    for (i=0; i<thisSG->dim; i++) fprintf (out," %g",currn->x[i]);
    fprintf (out,"\n");
    currn->index = ++nnode;
    currn = currn->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // then, all the radii
  currr = thisSG->radii->first;
  while (currr) {
    fprintf (out,"vr %g\n",currr->r);
    currr->index = ++nrad;
    currr = currr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }
*/

  // then, all the tangents
  currt = thisSG->tangents->first;
  while (currt) {
    fprintf (out,"vn");
    for (i=0; i<thisSG->dim; i++) fprintf (out," %g",currt->x[i]);
    fprintf (out,"\n");
    currt->index = ++ntan;
    currt = currt->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // then all of the segments
  curr = thisSG->first;
  while (curr) {
    // print the nodes, if they've not yet been printed
    for (i=0; i<2; i++) {
      if (curr->n[i]->flag) {
        fprintf (out,"v");
        for (j=0; j<thisSG->dim; j++) fprintf (out," %g",curr->n[i]->x[j]);
        fprintf (out,"\n");
        curr->n[i]->index = ++nnode;
        curr->n[i]->flag = false;
      }
    }
    // print the radius, if it hasn't been printed yet
    for (i=0; i<2; i++) {
      if (curr->r[i]) {
      if (curr->r[i]->flag) {
        fprintf (out,"vr %g\n",curr->r[i]->r);
        curr->r[i]->index = ++nrad;
        curr->r[i]->flag = false;
      }
      }
    }
    // then, print the segment
    fprintf (out,"s");
    for (i=0; i<2; i++) {
      if (curr->r[i]) {
        if (curr->t[i]) {
          fprintf (out," %ld/%ld/%ld",curr->n[i]->index,curr->r[i]->index,
                                      curr->t[i]->index);
        } else {
          fprintf (out," %ld/%ld",curr->n[i]->index,curr->r[i]->index);
        }
      } else {
        if (curr->t[i]) {
          fprintf (out," %ld//%ld",curr->n[i]->index,curr->t[i]->index);
        } else {
          fprintf (out," %ld",curr->n[i]->index);
        }
      }
    }
    fprintf (out,"\n");
    curr = curr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  fprintf(stderr,"%d elements\n",nlines);
  fflush(stderr);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Write a .vtk file (steal syntax from fbbconvert.c)
 */
int write_vtk (FILE *out, seg_group_ptr thisSG) {

  int i,j;
  int thisdim = 3;
  int nnode = 0;
  int nrad = 0;
  int ntan = 0;
  int nlines = 0;
  seg_ptr curr;

  fprintf (stderr,"Writing .vtk file");
  fflush (stderr);

  if (thisSG->dim < thisdim) thisdim = thisSG->dim;

  // for this first take, blindly dump all segments and nodes

  fprintf (out,"# vtk DataFile Version 2.0\n");
  fprintf (out,"%s\n",VERSION);
  fprintf (out,"ASCII\n");
  fprintf (out,"DATASET UNSTRUCTURED_GRID\n");
  // note: we assume that the node count is correct
  fprintf (out,"POINTS %d float\n",thisSG->nodes->num);
  fprintf(stderr,"."); fflush(stderr);

  // first, print the nodes
  nnode = write_vtk_nodes_and_set_indexes (out, thisSG->nodes, thisdim, 0);
  fprintf(stderr,"."); fflush(stderr);

/*
  // then, all the radii
  currr = thisSG->radii->first;
  while (currr) {
    fprintf (out,"vr %g\n",currr->r);
    currr->index = ++nrad;
    currr = currr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  // then, all the tangents
  currt = thisSG->tangents->first;
  while (currt) {
    fprintf (out,"vn");
    for (i=0; i<thisSG->dim; i++) fprintf (out," %g",currt->t[i]);
    fprintf (out,"\n");
    currt->index = ++ntan;
    currt = currt->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }
*/

  // then all of the segments
  fprintf (out,"CELLS %d %d\n",thisSG->num,thisSG->num*3);
  curr = thisSG->first;
  while (curr) {
    // then, print the segment
    fprintf (out,"2 %ld %ld\n",curr->n[0]->index,curr->n[1]->index);
    curr = curr->next;
  }
  fprintf(stderr,"."); fflush(stderr);

  fprintf (out,"CELL_TYPES %d\n",thisSG->num);
  for (i=0; i<thisSG->num; i++) fprintf (out,"3\n");
  fprintf(stderr,"."); fflush(stderr);


  // then, scalar values, namely, the radii
  fprintf (out,"POINT_DATA %d\n",thisSG->nodes->num);
  fprintf (out,"SCALARS radius float 1\n");
  fprintf (out,"LOOKUP_TABLE default\n");
  (void) write_vtk_nodes_radii (out, thisSG->nodes, thisSG->radius);


  fprintf(stderr,".\n"); fflush(stderr);

  //fprintf(stderr,"%d elements\n",nlines);
  //fflush(stderr);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Recursively index and print nodes to .vtk file
 *
 * Note: indexes start at 0 !
 */
int write_vtk_nodes_and_set_indexes (FILE *outfp, node_group_ptr thisNG,
  int dim, int sum) {

  int i;
  node_ptr curr;

  // first, see if we have this entry already
  if (thisNG->child[0]) {

    // there are children, call each of them
    sum = write_vtk_nodes_and_set_indexes (outfp, thisNG->child[0], dim, sum);
    sum = write_vtk_nodes_and_set_indexes (outfp, thisNG->child[1], dim, sum);

  } else {

    // there are no children, check all local entries
    curr = thisNG->first;
    while (curr) {
      fprintf (outfp,"%g",curr->x[0]);
      for (i=1; i<dim; i++) fprintf (outfp," %g",curr->x[i]);
      fprintf (outfp,"\n");
      curr->index = sum++;
      curr = curr->next;
    }
  }

  return (sum);
}


/*
 * Recursively index and print nodes to .vtk file
 *
 * Note: indexes start at 0 !
 */
int write_vtk_nodes_radii (FILE *outfp, node_group_ptr thisNG, double rad) {

  int i;
  double bigrad;
  node_ptr curr;

  // first, see if we have this entry already
  if (thisNG->child[0]) {

    // there are children, call each of them
    (void) write_vtk_nodes_radii (outfp, thisNG->child[0], rad);
    (void) write_vtk_nodes_radii (outfp, thisNG->child[1], rad);

  } else {

    // there are no children, check all local entries
    curr = thisNG->first;
    while (curr) {

      //fprintf (stderr, "node at %g %g %g\n",curr->x[0],curr->x[1],curr->x[2]);

      //for (i=0; i<curr->numconn0; i++) 
      //fprintf (stderr, "  0 conn %d has nodes %ld %ld\n",i,
      //  curr->conn0[i]->n[0]->index,curr->conn0[i]->n[1]->index);
      //for (i=0; i<curr->numconn1; i++) 
      //fprintf (stderr, "  1 conn %d has nodes %ld %ld\n",i,
      //  curr->conn1[i]->n[0]->index,curr->conn1[i]->n[1]->index);

      // choose the largest radius at this node
      bigrad = -1.;
      for (i=0; i<curr->numconn0; i++) 
        if (curr->conn0[i]->r[0])
          if (curr->conn0[i]->r[0]->r > bigrad)
            bigrad = curr->conn0[i]->r[0]->r;
      //fprintf (stderr, "  rad %g\n",bigrad);
      for (i=0; i<curr->numconn1; i++) 
        if (curr->conn1[i]->r[1])
          if (curr->conn1[i]->r[1]->r > bigrad)
            bigrad = curr->conn1[i]->r[1]->r;
      //fprintf (stderr, "  rad %g\n",bigrad);
      // print either the largest radius, or the default radius
      if (bigrad > -EPSILON) fprintf (outfp,"%g\n",bigrad);
      else fprintf (outfp,"%g\n",rad);

      curr = curr->next;
    }
  }

  return (0);
}


/*
 * Write a .seg file (syntax very much like .obj)
 */
int write_svg (FILE *out, seg_group_ptr thisSG, int argc, char **argv) {

  fprintf (stderr,"Writing .svg file...");
  fflush (stderr);
  static int fileCnt = 0;

  // set all flags to true (true==not been printed yet)
  //(void) set_radius_flags (thisSG->radii, true);
  //(void) set_node_flags (thisSG->nodes, true);

  // make a file name
  //char outfile[255];
  //sprintf(outfile,"out_%04d.svg",fileCnt);

  // call connected C++ routine for this work
  int retVal = write_svg_using_wxSVG(out, thisSG, argc, argv);

  fprintf(stderr,"%d elements\n",retVal);
  fflush(stderr);

  fileCnt++;

  return(retVal);
}


/*
 * Recursively set all flags
 */
int set_seg_flags (seg_group_ptr thisSG, const bool val) {

  seg_ptr curr;

  // there are no children, check all local entries
  curr = thisSG->first;
  while (curr) {
    curr->flag = val;
    curr = curr->next;
  }

  return (0);
}


/*
 * Recursively set all flags
 */
int set_node_flags (node_group_ptr thisNG, const bool val) {

  node_ptr curr;

  // first, see if we have this entry already
  if (thisNG->child[0]) {

    // there are children, call each of them
    (void) set_node_flags (thisNG->child[0], val);
    (void) set_node_flags (thisNG->child[1], val);

  } else {

    // there are no children, check all local entries
    curr = thisNG->first;
    while (curr) {
      curr->flag = val;
      curr = curr->next;
    }
  }

  return (0);
}


/*
 * Recursively set all flags
 */
int set_radius_flags (rad_group_ptr thisRG, const bool val) {

  rad_ptr curr;

  // first, see if we have this entry already
  if (thisRG->child[0]) {

    // there are children, call each of them
    (void) set_radius_flags (thisRG->child[0], val);
    (void) set_radius_flags (thisRG->child[1], val);

  } else {

    // there are no children, check all local entries
    curr = thisRG->first;
    while (curr) {
      curr->flag = val;
      curr = curr->next;
    }
  }

  return (0);
}


/*
 * Write a Radiance file (fat)
 */
int write_rad (FILE *out, seg_group_ptr thisSG) {

  int write_no_cones = false;
  int write_many_colors = false;
  int i;
  int count = 0;
  int nlines = 0;
  node_ptr currn;
  seg_ptr curr;
  int thisdim = thisSG->dim;
  double thisrad = -1.;
  char color[5];

  fprintf (stderr,"Writing .rad file");
  fflush (stderr);
  strcpy(color,"def");

  // make sure we dump 3-dimensional data
  thisdim = thisSG->dim;
  if (thisdim > 3) thisdim = 3;

  // for this first take, blindly dump all segments and nodes
  fprintf (out,"# %s\n",VERSION);

  // first, the nodes
  (void) write_radiance_nodes (out, thisSG->nodes, thisdim, thisSG->radius,
                               &count, &nlines);

  // then, all of the segments
  count = 0;
  curr = thisSG->first;
  while (curr) {

    if (write_many_colors) {
      //sprintf(color,"def%1d",(int)(10.*rand()/(double)(RAND_MAX)));
      sprintf(color,"def%1d",(int)(-0.5*log(rand()/(double)(RAND_MAX))));
    }

    // do we write a cylinder or a cone?
    if (curr->r[0] && curr->r[1]) {
      // if we have radii for both nodes...
      if (write_no_cones) {
        // find smallest radius
        if (curr->r[0]->r > curr->r[1]->r) thisrad = curr->r[0]->r;
        else thisrad = curr->r[1]->r;
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",thisrad);
      } else if (fabs(curr->r[0]->r - curr->r[1]->r) < EPSILON) {
        // and they're the same radius, write a cylinder
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",curr->r[0]->r);
      } else {
        // different radius, write a cone
        fprintf (out,"%s cone c%ld 0 0 8",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e %12.8e\n",curr->r[0]->r,curr->r[1]->r);
      }
    } else if (curr->r[0]) {
      // if we have radii for only one node...
      if (write_no_cones) {
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",curr->r[0]->r);
      } else if (fabs(curr->r[0]->r - thisSG->radius) < EPSILON) {
        // same radius, write a cylinder
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",thisSG->radius);
      } else {
        // different radius, write a cone
        fprintf (out,"%s cone c%ld 0 0 8",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e %12.8e\n",curr->r[0]->r,thisSG->radius);
      }
    } else if (curr->r[1]) {
      // if we have radii for only one node...
      if (write_no_cones) {
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",curr->r[1]->r);
      } else if (fabs(curr->r[1]->r - thisSG->radius) < EPSILON) {
        // same radius, write a cylinder
        fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",thisSG->radius);
      } else {
        // different radius, write a cone
        fprintf (out,"%s cone c%ld 0 0 8",color,++count);
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e %12.8e\n",thisSG->radius,curr->r[0]->r);
      }
    } else {
      // if we have no radii for either node...
      // use the global radius and write a cylinder
      fprintf (out,"%s cylinder c%ld 0 0 7",color,++count);
      for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
      for (i=thisdim; i<3; i++) fprintf (out," 0.");
      for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
      for (i=thisdim; i<3; i++) fprintf (out," 0.");
      fprintf (out," %12.8e\n",thisSG->radius);
    }

    curr = curr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  fprintf(stderr,"%d elements\n",nlines);
  fflush(stderr);

  // if things bomb, send back a nonzero
  return (0);
}


/*
 * Recursively write the Radiance spheres
 */
int write_radiance_nodes (FILE *out, node_group_ptr thisNG, int dim,
  double rad, int *count, int *nlines) {

  int write_many_colors = false;
  int i;
  double bigrad;
  node_ptr curr;
  char color[5];

  strcpy(color,"def");

  if (thisNG->child[0]) {

    (void) write_radiance_nodes (out, thisNG->child[0], dim, rad, count, nlines);
    (void) write_radiance_nodes (out, thisNG->child[1], dim, rad, count, nlines);

  } else {

    curr = thisNG->first;
    while (curr) {
      // choose the largest radius at this node
      bigrad = -1.;
      for (i=0; i<curr->numconn0; i++) 
        if (curr->conn0[i]->r[0])
          if (curr->conn0[i]->r[0]->r > bigrad)
            bigrad = curr->conn0[i]->r[0]->r;
      for (i=0; i<curr->numconn1; i++) 
        if (curr->conn1[i]->r[1])
          if (curr->conn1[i]->r[1]->r > bigrad)
            bigrad = curr->conn1[i]->r[1]->r;
      // print either the largest radius, or the default radius
      if (bigrad < -EPSILON) bigrad = rad;
      // but only print a sphere if it's large enough!
      if (bigrad > EPSILON) {
        if (write_many_colors) {
          sprintf(color,"def%1d",(int)(-0.5*log(rand()/(double)(RAND_MAX))));
        }
        fprintf (out,"%s sphere s%ld 0 0 4",color,++(*count));
        for (i=0; i<dim; i++) fprintf (out," %12.8e",curr->x[i]);
        for (i=dim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e\n",bigrad);
      }

      curr = curr->next;

      if (++(*nlines)%DOTPER == 1) {
        fprintf(stderr,".");
        fflush(stderr);
      }
    }
  }

  return (0);
}


/*
 * Calculate and write the appropriate dots for a given structure
 */
int write_dots (FILE *out, seg_group_ptr thisSG, double res) {

  int i;
  int write_vtk_headers = false;
  double *thresh;
  node_group_ptr dots;

  // make a structure for the dots to print
  dots = (NODE_GROUP*) malloc (sizeof(NODE_GROUP));
  dots->num = 0;
  dots->first = NULL;
  dots->parent = NULL;
  dots->child[0] = NULL;
  dots->child[1] = NULL;
  dots->axis = -1;
  dots->min = NULL;
  dots->max = NULL;

  // set the dot proximity threshold
  thresh = (double*) malloc (thisSG->dim*sizeof(double));
  for (i=0; i<thisSG->dim; i++) thresh[i] = res;

  // begin filling the dots structure
  (void) write_dots_1 (thisSG->nodes, dots, thisSG->dim, thresh);
  fprintf(stderr,"  made %d dots from nodes\n",dots->num);

  // now, loop over all segments, adding any dots that will fit
  (void) write_dots_2 (thisSG, dots, thisSG->dim, thresh);
  fprintf(stderr,"  now %d dots after including segments\n",dots->num);

  // finally, write the data
  if (write_vtk_headers) {
    fprintf (out,"# vtk DataFile Version 2.0\n");
    fprintf (out,"%s\n",VERSION);
    fprintf (out,"ASCII\n");
    fprintf (out,"DATASET UNSTRUCTURED_GRID\n");
    // note: we assume that the node count is correct
    fprintf (out,"POINTS %d float\n",dots->num);
    // finally, loop through the list of dots, writing what we've found
  }

  // write all of the dot locations
  (void) write_dots_3 (out, dots, thisSG->dim);

  // then the rest of the data
  if (write_vtk_headers) {
    fprintf (out,"CELLS %d %d\n",dots->num,dots->num*2);
    for (i=0; i<dots->num; i++) fprintf (out,"1 %d\n",i);
    fprintf (out,"CELL_TYPES %d\n",dots->num);
    for (i=0; i<dots->num; i++) fprintf (out,"1\n");
  }

  // free the tree structure and all of the dots
  free_nodes (dots);

  return (0);
}


/*
 * Recursively include nodes as dots
 */
int write_dots_1 (node_group_ptr thisNG, node_group_ptr dots,
  int dim, double *thresh) {

  int i;
  double dist;
  node_ptr curr, close_dot;

  if (thisNG->child[0] != NULL) {

    (void) write_dots_1 (thisNG->child[0], dots, dim, thresh);
    (void) write_dots_1 (thisNG->child[1], dots, dim, thresh);

  } else {

    curr = thisNG->first;
    while (curr != NULL) {

      // set the threshold distance
      dist = thresh[0];

      // check to see if this node is too close to any other node
      close_dot = find_closest_node (dots, dim, curr->x, NULL, &dist);
      //close_dot = NULL;
      //close_dot = curr;

      // if it is not too close, add it to the dots list
      if (close_dot == NULL) {
        (void) add_node (dots, dim, curr->x, 0);
        //fprintf(stderr,"  added dot, closest node was %g\n",dist);
      //} else {
        //fprintf(stderr,"    no dot, closest node was %g\n",dist);
      }

      curr = curr->next;
    }
  }

  return (0);
}


/*
 * Recursively include segments as dots
 */
int write_dots_2 (seg_group_ptr thisSG, node_group_ptr dots,
  int dim, double *thresh) {

  int i,isubseg,nsubseg;
  double dist,len;
  double *testloc, *tangent;
  seg_ptr curr;
  node_ptr close_dot;

  //if (thisSG->child[0] != NULL) {

    //(void) write_dots_2 (thisSG->child[0], dots, dim, thresh);
    //(void) write_dots_2 (thisSG->child[1], dots, dim, thresh);

  //} else {

    testloc = (double*) malloc (dim*sizeof(double));
    tangent = (double*) malloc (dim*sizeof(double));

    curr = thisSG->first;
    while (curr != NULL) {

      // how long is this segment?
      len = seg_length (curr, dim);
      nsubseg = (int)(len/thresh[0]);

      // find the straight tangential vector along this edge
      for (i=0; i<dim; i++) tangent[i] = curr->n[1]->x[i]-curr->n[0]->x[i];

      // march along the segment, checking locations for particles
      for (isubseg=1; isubseg<nsubseg; isubseg++) {

        // set the test location
        for (i=0; i<dim; i++) testloc[i] = curr->n[0]->x[i] +
          (double)isubseg * tangent[i] / (double)nsubseg;

        // check to see if this dot is too close to any other dot
        dist = thresh[0];
        close_dot = find_closest_node (dots, dim, testloc, NULL, &dist);

        // if it is not too close, add it to the dots list
        if (close_dot == NULL) (void) add_node (dots, dim, testloc, 0);
      }

      // if the segment is between 1 and 2 thresholds long, try two points
      //if (nsubseg == 1) {
      if (false) {

        // set the test location
        for (i=0; i<dim; i++) testloc[i] = curr->n[0]->x[i] +
          (thresh[0]*1.01) * tangent[i] / len;

        // check to see if this dot is too close to any other dot
        dist = thresh[0];
        close_dot = find_closest_node (dots, dim, testloc, NULL, &dist);

        // if it is not too close, add it to the dots list
        if (close_dot == NULL) (void) add_node (dots, dim, testloc, 0);

        // now, try the dot thresh distance from the other node
        for (i=0; i<dim; i++) testloc[i] = curr->n[1]->x[i] -
          (thresh[0]*1.01) * tangent[i] / len;
        dist = thresh[0];
        close_dot = find_closest_node (dots, dim, testloc, NULL, &dist);
        if (close_dot == NULL) (void) add_node (dots, dim, testloc, 0);

      }

      curr = curr->next;
    }
  //}

  return (0);
}


/*
 * Recursively write dots as ASCII file
 */
int write_dots_3 (FILE *out, node_group_ptr thisNG, int dim) {

  int i;
  node_ptr curr;

  if (thisNG->child[0] != NULL) {

    (void) write_dots_3 (out, thisNG->child[0], dim);
    (void) write_dots_3 (out, thisNG->child[1], dim);

  } else {

    curr = thisNG->first;
    while (curr != NULL) {

      // dump the data
      fprintf(out,"%g",curr->x[0]);
      for (i=1; i<dim; i++) fprintf(out," %g",curr->x[i]);
      fprintf(out,"\n");

      curr = curr->next;
    }
  }

  return (0);
}


/*
 * March along paths and write triangle mesh
 */
int write_obj (FILE *out, seg_group_ptr segs, const double mindx) {

  int i;
  int count = 0;
  int nlines = 0;
  seg_ptr curr;

  fprintf (stderr,"Writing .obj file\n");
  fflush (stderr);

  // first, identify separate disconnected strands
  const int nstrands = identify_separate_strands (segs);
  fprintf (stderr,"  found %d separate strands\n", nstrands);

  // set flags to indicate which segments are done
  (void) set_seg_flags (segs, false);
  // maybe set the node flags, too?
  (void) set_node_flags (segs->nodes, false);

  // make sure we dump 3-dimensional data (add 0.0 for z dim if 2D, ignore dims higher than 3)
  int thisdim = segs->dim;
  if (thisdim > 3) thisdim = 3;

  // set tangent vectors on all segments
  set_seg_tangents (segs, thisdim);

  // the current and previous basis vectors
  // b1 = tangent along strand
  // b2 = first perpendicular axis
  // b3 = second perpendicular axis = b1 x b2
  double thisb1[3], thisb2[3], thisb3[3];
  double lastb1[3], lastb2[3], lastb3[3];

  // these need to be persistent across strands
  int lastidx = 1;
  int lastcnt = 0;
  int curridx = -1;
  int currcnt = -1;

  // for this first take, blindly dump all segments and nodes
  fprintf (out,"# %s\n",VERSION);

  // write one strand at a time
  for (int istrand=0; istrand<nstrands; ++istrand) {

    // find one end of this strand
    node_ptr currn = find_end_node_of_strand(segs, istrand);
    node_ptr lastn = NULL;
    node_ptr nextn = NULL;
    seg_ptr currseg = NULL;
    seg_ptr nextseg = NULL;

    while (currn) {
      fprintf (stderr,"node at %g %g %g\n", currn->x[0], currn->x[1], currn->x[2]);

      // find next node and next segment
      if (currn->numconn0 > 0) {
        if (! currn->conn0[0]->n[1]->flag) {
          nextseg = currn->conn0[0];
          nextn = currn->conn0[0]->n[1];
          fprintf (stderr,"    next node is %d along segment %d\n", nextn->index, nextseg->index);
        } else if (currn->numconn1 > 0) {
          if (! currn->conn1[0]->n[0]->flag) {
            nextseg = currn->conn1[0];
            nextn = currn->conn1[0]->n[0];
            fprintf (stderr,"    next node is %d along segment %d\n", nextn->index, nextseg->index);
          }
        } else {
          nextn = NULL;
          nextseg = NULL;
        }
      } else {
        nextn = NULL;
        nextseg = NULL;
      }

      // ensure we have a radius for each end
      double r0 = segs->radius;
      // find which radius to use!
      //if (curr->r[0]) r0 = curr->r[0]->r;
      // how many nodes on each ring?
      //int n0 = (int)(0.5 + (6.28318531 * r0) / mindx);
      //if (n0 < 3) n0 = 3;
      //if (r0 == 0.0) n0 = 1;
      int n0 = 13;

      // set this node's basis vectors - first the tangent
      seg_ptr testseg = NULL;
      if (currseg) testseg = currseg;
      if (nextseg) testseg = nextseg;
      tan_ptr tptr = NULL;
      if (testseg->n[0] == currn) {
        tptr = testseg->t[0];
      } else {
        tptr = testseg->t[1];
      }
      for (i=0; i<3; i++) thisb1[i] = tptr->x[i];
      vec_normalize_in_place(thisb1, 3);
      fprintf (stderr,"    tangent is %g %g %g\n", thisb1[0], thisb1[1], thisb1[2]);

      // now set perpendicular basis vectors!
      if (lastn == NULL) {
        // if this is the first node in the strand, make one up!
        if (abs(thisb1[0]) > 0.9) {
          thisb2[0] = 0.0; thisb2[1] = 1.0; thisb2[2] = 0.0;
        } else {
          thisb2[0] = 1.0; thisb2[1] = 0.0; thisb2[2] = 0.0;
        }
      } else {
        // this is not the first node, so rotate the two other basis vectors
        // first, find the axis and angle of rotation from node lastn to currn
        const double angle = acos(vec_dot(lastb1, thisb1, 3));
        fprintf (stderr,"    angle is %g\n", angle);
        // acos always returns 0..pi
        if (angle < 0.01) {
          // if the angle is minuscule, just copy the 2nd basis vector
          for (i=0; i<3; i++) thisb2[i] = lastb2[i];
        } else {
          // if the angle is significant, must rotate b2 around the rotational axis
          double axis[3];
          vec_cross (lastb1, thisb1, axis);
          vec_normalize_in_place(axis, 3);
          fprintf (stderr,"    rotation is %g %g %g\n", axis[0], axis[1], axis[2]);
          // use a rotation matrix from axis-angle
          const double ca = 1.0 - cos(angle);
          const double omc = 1.0 - cos(angle);
          const double sa = sin(angle);
          double newb2[3];
          newb2[0] = lastb2[0] * (ca+axis[0]*axis[0]*omc) +
                     lastb2[1] * (axis[0]*axis[1]*omc-axis[2]*sa) +
                     lastb2[2] * (axis[0]*axis[2]*omc+axis[1]*sa);
          newb2[1] = lastb2[0] * (axis[1]*axis[0]*omc+axis[2]*sa) +
                     lastb2[1] * (ca+axis[1]*axis[1]*omc) +
                     lastb2[2] * (axis[1]*axis[2]*omc-axis[0]*sa);
          newb2[2] = lastb2[0] * (axis[2]*axis[0]*omc-axis[1]*sa) +
                     lastb2[1] * (axis[2]*axis[1]*omc+axis[0]*sa) +
                     lastb2[2] * (ca+axis[2]*axis[2]*omc);
        }
      }
      // make the first vector perpendicular to the second one
      vec_make_perpendicular(thisb2, thisb1, 3);
      vec_normalize_in_place(thisb2, 3);
      // calc b3 from b1, b2
      vec_cross (thisb1, thisb2, thisb3);
      fprintf (stderr,"    basis2 is %g %g %g\n", thisb2[0], thisb2[1], thisb2[2]);
      fprintf (stderr,"    basis3 is %g %g %g\n", thisb3[0], thisb3[1], thisb3[2]);

      // if this is the first node AND the radius is not 0, get ready for an endcap
      if (lastn == NULL && n0 > 1) {
        fprintf (out,"v");
        for (i=0; i<thisdim; i++) fprintf (out," %g", currn->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out,"\n");
        lastidx = lastidx + lastcnt;
        lastcnt = 1;
      }

      // write the ring of output nodes
      for (int j=0; j<n0; ++j) {
        const double theta = 2. * M_PI * (j+0.5) / (double)n0;
        const double st = r0 * sin(theta);
        const double ct = r0 * cos(theta);
        fprintf (out,"v");
        for (i=0; i<thisdim; i++) fprintf (out," %g", currn->x[i] + ct*thisb2[i] + st*thisb3[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out,"\n");
        //fprintf (out," # node %d\n", currn->index);
      }
      // index of the first vertex (1-indexed)
      curridx = lastidx + lastcnt;
      currcnt = n0;
      fprintf (stderr,"    last nodes are %d to %d\n", lastidx, curridx-1);
      fprintf (stderr,"    curr nodes are %d to %d\n", curridx, curridx+currcnt-1);

      // write any triangles necessary
      if (currseg) {
        for (int j=0; j<n0-1; ++j) {
          fprintf (out,"f %d %d %d\n", lastidx+j, lastidx+j+1, curridx+j);
          fprintf (out,"f %d %d %d\n", lastidx+j+1, curridx+j+1, curridx+j);
        }
        // the last pair needs special treatment
        fprintf (out,"f %d %d %d\n", lastidx+lastcnt-1, lastidx, curridx+currcnt-1);
        fprintf (out,"f %d %d %d\n", lastidx, curridx, curridx+currcnt-1);
      } else if (lastn == NULL && n0 > 1) {
        fprintf (stderr,"    making endcap around node %d\n", lastidx);
        // this is the starting endcap
        for (int j=0; j<n0-1; ++j) {
          fprintf (out,"f %d %d %d\n", curridx+j, lastidx, curridx+j+1);
        }
        // the last one needs special treatment
        fprintf (out,"f %d %d %d\n", curridx+currcnt-1, lastidx, curridx);
      }

      // make the ending endcap
      if (nextseg == NULL && n0 > 1) {
        fprintf (stderr,"    making endcap around node %d\n", curridx+currcnt);
        // first the central vertex
        fprintf (out,"v");
        for (i=0; i<thisdim; i++) fprintf (out," %g", currn->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out,"\n");
        // then the triangles
        for (int j=0; j<n0-1; ++j) {
          fprintf (out,"f %d %d %d\n", curridx+j, curridx+j+1, curridx+currcnt);
        }
        // the last one needs special treatment
        fprintf (out,"f %d %d %d\n", curridx+currcnt-1, curridx, curridx+currcnt);
        // update currcnt for the next ring
        currcnt += 1;
      }

      // mark this node as done
      currn->flag = true;

      //for (int i=0; i<currn->numconn0; ++i) {
      //  seg_ptr thisseg = currn->conn0[i];
      //  fprintf (stderr,"    conn0 seg %d has nodes %d %d\n", thisseg->index, thisseg->n[0]->index, thisseg->n[1]->index);
      //}
      //for (int i=0; i<currn->numconn1; ++i) {
      //  seg_ptr thisseg = currn->conn1[i];
      //  fprintf (stderr,"    conn1 seg %d has nodes %d %d\n", thisseg->index, thisseg->n[0]->index, thisseg->n[1]->index);
      //}

      // copy current to last
      lastn = currn;
      currseg = nextseg;
      currn = nextn;
      lastidx = curridx;
      lastcnt = currcnt;
      for (i=0; i<3; i++) lastb1[i] = thisb1[i];
      for (i=0; i<3; i++) lastb2[i] = thisb2[i];
      for (i=0; i<3; i++) lastb3[i] = thisb3[i];
    }

  }
  return 1;

  // then, all of the segments
  count = 0;
  curr = segs->first;
  while (curr) {

    // ensure we have a radius for each end
    double r0 = segs->radius;
    if (curr->r[0]) r0 = curr->r[0]->r;
    double r1 = segs->radius;
    if (curr->r[1]) r1 = curr->r[1]->r;

    // how many nodes on each ring?
    int n0 = (int)(0.5 + (6.28318531 * r0) / mindx);
    if (n0 < 3) n0 = 3;
    if (r0 == 0.0) n0 = 1;
    int n1 = (int)(0.5 + (6.28318531 * r1) / mindx);
    if (n1 < 3) n1 = 3;
    if (r1 == 0.0) n1 = 1;

    // do we write a cylinder or a cone?
    if (curr->r[0] && curr->r[1]) {
      // if we have radii for both nodes...
      {
        // different radius, write a cone
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
        for (i=thisdim; i<3; i++) fprintf (out," 0.");
        fprintf (out," %12.8e %12.8e\n",curr->r[0]->r,curr->r[1]->r);
      }
    } else {
      // if we have no radii for either node...
      // use the global radius and write a cylinder
      for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[0]->x[i]);
      for (i=thisdim; i<3; i++) fprintf (out," 0.");
      for (i=0; i<thisdim; i++) fprintf (out," %12.8e",curr->n[1]->x[i]);
      for (i=thisdim; i<3; i++) fprintf (out," 0.");
      fprintf (out," %12.8e\n",segs->radius);
    }

    // don't just use whichever node is next in the list!
    // must use the adjacent segment instead!
    curr = curr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  fprintf(stderr,"%d elements\n",nlines);
  fflush(stderr);

  // if things bomb, send back a nonzero
}


/*
 * Create a Gaussian perturbation in as many dimensions as required,
 * mean is 0, std dev is scale
 */
void perturb_gaussian (double *dx, double scale, int dim) {

  int d;
  double t1,t2;

  // do two at a time
  for (d=0; d<dim-1; d+=2) {
    t1 = (double)(rand())/(double)(RAND_MAX);
    t2 = (double)(rand())/(double)(RAND_MAX);
    dx[d] = scale*sqrt(-2.*log(t1))*cos(2.*M_PI*t2);
    dx[d+1] = scale*sqrt(-2.*log(t1))*sin(2.*M_PI*t2);
  }

  // do the last one (for odd-numbered dimensions)
  if (d != dim) {
    t1 = (double)(rand())/(double)(RAND_MAX);
    t2 = (double)(rand())/(double)(RAND_MAX);
    dx[dim-1] = scale*sqrt(-2.*log(t1))*cos(2.*M_PI*t2);
  }

  return;
}


/* distance */
double seg_length (seg_ptr thisSP, int dim) {
  int i;
  double vds = 0.;
  for (i=0; i<dim; i++) vds += pow(thisSP->n[0]->x[i]-thisSP->n[1]->x[i],2);
  return (sqrt(vds));
}


/* distance squared */
double seg_length_sqrd (seg_ptr thisSP, int dim) {
  int i;
  double vds = 0.;
  for (i=0; i<dim; i++) vds += pow(thisSP->n[0]->x[i]-thisSP->n[1]->x[i],2);
  return (vds);
}


/* distance squared */
double vec_dist_sqrd (double *x, double *y, int dim) {
  int i;
  double vds = 0.;
  for (i=0; i<dim; i++) vds += pow(x[i]-y[i],2);
  return (vds);
}


/* inverse square root length */
double vec_inv_sqrt (double *x, int dim) {
  int i;
  double vds = 0.;
  for (i=0; i<dim; i++) vds += pow(x[i],2);
  return (1./sqrt(vds));
}


/* dot product */
double vec_dot (double const * const x, double const * const y, const int dim) {
  double dot = 0.;
  for (int i=0; i<dim; i++) dot += x[i]*y[i];
  return dot;
}


/* make the first vector perpendicular to the second one */
void vec_make_perpendicular (double * const x, double const * const dir, const int dim) {
  // first, find the length of the direction vector
  double dirlen = 0.;
  for (int i=0; i<dim; i++) dirlen += dir[i]*dir[i];
  dirlen = 1./sqrt(dirlen);

  // now, find the dot product of the normalized direction vector and the test vector
  double dotp = 0.;
  for (int i=0; i<dim; i++) dotp += x[i]*dir[i]*dirlen;

  // finally, subtract the part of x that is along dir
  for (int i=0; i<dim; i++) x[i] -= dotp*dir[i]*dirlen;

  return;
}

/* normalize a vector to unit length */
void vec_normalize (double *x, double *xnorm, const int dim) {
  int i;
  double dotp = 0.;

  for (i=0; i<dim; i++) dotp += x[i]*x[i];
  dotp = 1./sqrt(dotp);
  for (i=0; i<dim; i++) xnorm[i] = dotp*x[i];

  return;
}

void vec_normalize_in_place (double *x, const int dim) {
  int i;
  double dotp = 0.;

  for (i=0; i<dim; i++) dotp += x[i]*x[i];
  dotp = 1./sqrt(dotp);
  for (i=0; i<dim; i++) x[i] *= dotp;

  return;
}

/* 3d vector cross product */
void vec_cross (double *x, double *y, double *prod) {
  prod[0] = x[1]*y[2] - y[1]*x[2];
  prod[1] = x[2]*y[0] - y[2]*x[0];
  prod[2] = x[0]*y[1] - y[0]*x[1];
  return;
}


/*
 * Recursively free memory associated with nodes
 */
int free_nodes (node_group_ptr thisNG) {

  node_ptr curr, next;

  if (thisNG->child[0] != NULL) {

    (void) free_nodes (thisNG->child[0]);
    (void) free_nodes (thisNG->child[1]);

  } else {

    curr = thisNG->first;
    while (curr != NULL) {
      free (curr->x);
      free (curr->conn0);
      free (curr->conn1);
      next = curr->next;
      free (curr);
      curr = next;
    }
  }

  free (thisNG->min);
  free (thisNG->max);
  free (thisNG);

  return (0);
}


/*
 * Recursively free memory associated with tangents
 */
int free_tangents (tan_group_ptr thisTG) {

  tan_ptr curr, next;

  if (thisTG->child[0] != NULL) {

    (void) free_tangents (thisTG->child[0]);
    (void) free_tangents (thisTG->child[1]);

  } else {

    curr = thisTG->first;
    while (curr != NULL) {
      free (curr->x);
      next = curr->next;
      free (curr);
      curr = next;
    }
  }

  free (thisTG->min);
  free (thisTG->max);
  free (thisTG);

  return (0);
}


/*
 * Recursively free memory associated with radii
 */
int free_radii (rad_group_ptr thisRG) {

  rad_ptr curr, next;

  if (thisRG->child[0] != NULL) {

    (void) free_radii (thisRG->child[0]);
    (void) free_radii (thisRG->child[1]);

  } else {

    curr = thisRG->first;
    while (curr != NULL) {
      next = curr->next;
      free (curr);
      curr = next;
    }
  }

  free (thisRG);

  return (0);
}


/*
 * Free memory associated with segments (everything, basically)
 */
int free_segs (seg_group_ptr thisSG) {

  seg_ptr curr, next;

  // first, remove radii, tangents, and nodes
  (void) free_radii (thisSG->radii);
  (void) free_tangents (thisSG->tangents);
  (void) free_nodes (thisSG->nodes);

  // then, delete all segments
  curr = thisSG->first;
  while (curr != NULL) {
    next = curr->next;
    free (curr);
    curr = next;
  }

  // finally, free the structure itself
  free (thisSG);

  return (0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAXSTR],int status) {

  static char **cpp, *help_message[] =
  {
     "where [options] are one or more of the following:",
     " ",
     "   -seg        write a seg file (like wavefront obj) of the segments (default)",
     " ",
     "   -svg        write a scalable vector graphic image (SVG)",
     " ",
     "   -rad        write a Radiance-readable file of the cylinders",
     " ",
     "   -vtk        write a VTK-format file for paraview, etc.",
     " ",
     "   -png res    write a PNG image file with given maximum resolution",
     " ",
     "   -obj dx     write a wavefront-format triangle mesh with minimum edge length dx",
     " ",
     "   -bob dx     write a 3D brick-of-bytes with cell size dx",
     " ",
     "   -coarsen [l]  coarsen structure so that no element is shorter than",
     "               length l; default l=1",
     " ",
     "   -refine [l]  refine structure so that no element is longer than",
     "               length l; default l=1",
     " ",
     "   -srefine [l]  refine structure so that no element is longer than",
     "               length l; use splines to smooth curves; default l=1",
     " ",
     "   -roughen [s]  roughen structure by perturbing nodes, s is scale",
     "               factor; default s=1",
     " ",
     "   -prune [l]  prune any segments within l nodes from a tip;",
     "               default l=1",
     " ",
     "   -treeradius [s [dim [val]]]",
     "               set all radii to mimic woody plant growth,",
     "               relative strength is s * 10^6 * E / (rho * g * l),",
     "               next two entries define which nodes are roots:",
     "               '0 5.0' means whichever node is closest to the x=5 plane,",
     "               '2 0' means whichever node is closest to the z=0 plane;",
     "               default= 1.0 0 0.0",
     " ",
     "   -translate [x [y [z]]]",
     "               translate structure by vector x,y,z; default= 0,0,0",
     " ",
     "   -scale [x [y [z]]]",
     "               scale structure by magnitudes x,y,z; default= 1,1,1",
     "               mirroring can be done using negative values,",
     "               if fewer than (dimension) values are given, the",
     "               last value will be used for the rest of the dimensions",
     " ",
     "   -rscale s   scale all radii by the given factor",
     "   -rscale smin smax   scale all radii linearly to the given range",
     " ",
     "   -gr [rad]   sets global radius to rad, overwrites old radius; default=1",
     " ",
     "   -info       dump x,y,z,r min/max and other information",
     " ",
     "   -split      split into two .rad files along the longest dimension",
     " ",
     "   -zeroindexed  indicates that segment indices in input file start at",
     "               0 instead of 1",
     " ",
     "   -version    returns version information",
     " ",
     "   -help       returns this help information",
     " ",
     "The input file should be in seg, obj, rad, bob, or png format",
     " ",
     "Operations are done in the order that they appear on the command-",
     "line, so make sure to list your input file first!",
     " ",
     "Options may be abbreviated to an unambiguous length (duh).",
     "Output is always to stdout, so make sure to redirect (> or |)",
     NULL
  };

  fprintf(stderr,"Usage:\n\n  %s [infile | options]\n\n", progname);
  for (cpp = help_message; *cpp; cpp++)
    fprintf(stderr, "%s\n", *cpp);
    fflush(stderr);
  exit(status);
  return(0);
}

