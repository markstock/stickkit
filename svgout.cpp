/*
 * stickkit.c
 * svgout.cpp
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


// neccessary defines
#include "stickkit.h"
#include <svg.h>		// /usr/include/wxSVG/svg.h
#include <stdio.h>	// C-style file pointers


// prototypes
// need the "extern C" to tell the linker that C will be calling this!
extern "C" int write_svg_using_wxSVG (FILE*, seg_group_ptr, int, char**);


/*
 * Write a .seg file (syntax very much like .obj)
 *
 * Ideas from http://wxsvg.sourceforge.net/using.html
 * And docs at http://www.sourcecodebrowser.com/wxsvg/1.0b11/_c_s_s_style_declaration_8h.html#ad9b30a3628d9508eedb6e525fb38760bacffaa6c91cbe7b722fd13a5c5df0d069
 */
int write_svg_using_wxSVG (FILE* out, seg_group_ptr thisSG, int argc, char **argv) {

  int i,j;
  int nnode = 0;
  int nrad = 0;
  int ntan = 0;
  int nlines = 0;
  node_ptr currn;
  rad_ptr currr;
  tan_ptr currt;
  seg_ptr curr;
  wxSVGDocument* svgDoc = new wxSVGDocument;
  wxSVGLineElement* line = NULL;

  //wxSVGPoint* point = NULL;
  //wxSVGPolylineElement* pline = NULL;


  // find the maximum bounds of the object
  double maxBounds[thisSG->dim];
  double minBounds[thisSG->dim];
  for (i=0; i<thisSG->dim; i++) {
    maxBounds[i] = -9.e+9;
    minBounds[i] = 9.e+9;
  }
  double maxRad = -9.e+9;
  double minRad = 9.e+9;

  /*
  currn = thisSG->nodes->first;
  while (currn) {
    fprintf(stderr,"\nnode is at %g %g",currn->x[0],currn->x[1]);
    for (i=0; i<thisSG->dim; i++) {
      if (currn->x[i] > maxBounds[i]) maxBounds[i] = currn->x[i];
      if (currn->x[i] < minBounds[i]) minBounds[i] = currn->x[i];
    }
    currn = currn->next;
  }
  */

  curr = thisSG->first;
  while (curr) {
    // print the nodes, if they've not yet been printed
    for (i=0; i<2; i++) {
      //fprintf(stderr,"\nNode is at %g %g",curr->n[i]->x[0],curr->n[i]->x[1]);
      for (j=0; j<thisSG->dim; j++) {
        if (curr->n[i]->x[j] > maxBounds[j]) maxBounds[j] = curr->n[i]->x[j];
        if (curr->n[i]->x[j] < minBounds[j]) minBounds[j] = curr->n[i]->x[j];
      }
      if (curr->r[i]) {
        if (curr->r[i]->r > maxRad) maxRad = curr->r[i]->r;
        if (curr->r[i]->r < minRad) minRad = curr->r[i]->r;
      }
    }
    curr = curr->next;
  }

  /*
  // then, all the radii
  currr = thisSG->radii->first;
  while (currr) {
    if (currr->r > maxRad) maxRad = currr->r;
    if (currr->r < minRad) minRad = currr->r;
    currr = currr->next;
  }
  */

  // scale things up/down
  double scale = maxBounds[0]-minBounds[0];
  if (maxBounds[1]-minBounds[1] > scale) scale = maxBounds[1]-minBounds[1];
  scale = 1000.0 / scale;

  // Begin the document
  double w=1000,h=1000;
  wxSVGSVGElement* svgElement = new wxSVGSVGElement;
  //fprintf(stderr,"\nx bounds are %g %g\n",minBounds[0],maxBounds[0]);
  //fprintf(stderr,"\ny bounds are %g %g\n",minBounds[1],maxBounds[1]);
  //fprintf(stderr,"\nradius bounds are %g %g\n",minRad,maxRad);
  svgElement->SetWidth(w);
  svgElement->SetHeight(h);
  svgDoc->AppendChild(svgElement);


  // Create lines for all of the segments
  curr = thisSG->first;
  while (curr) {

    line = new wxSVGLineElement;
    line->SetStroke(wxSVGPaint(0,0,0));
    line->SetStrokeLinecap(wxCSS_VALUE_ROUND);

    // get the mean radius of this segment
    double radSum = 9.9e+9;
    if (curr->r[0]) {
      if (curr->r[0]->r < radSum) radSum = curr->r[0]->r;
    }
    if (curr->r[1]) {
      if (curr->r[1]->r < radSum) radSum = curr->r[1]->r;
    }
    if (radSum > 9.e+9) radSum = thisSG->radius;
    line->SetStrokeWidth(scale*radSum);

    // then, print the segment
    line->SetX1(scale*(curr->n[0]->x[0] - minBounds[0]));
    line->SetY1(scale*(curr->n[0]->x[1] - minBounds[1]));
    line->SetX2(scale*(curr->n[1]->x[0] - minBounds[0]));
    line->SetY2(scale*(curr->n[1]->x[1] - minBounds[1]));

    svgElement->AppendChild(line);

    curr = curr->next;

    if (++nlines%DOTPER == 1) {
      fprintf(stderr,".");
      fflush(stderr);
    }
  }

  //std::string outfile = outFileName;
  //svgDoc->Save(outfile);
  // make a file name
  //char outfile[255];
  //sprintf(outfile,"out_%04d.svg",fileCnt);

  //wxFile f (wxFile::wxFile::fd_stdout);
  //wxOuputStream fOutputStream (f);
  //svgDoc->Save(fOutputStream);
  //svgDoc->Save(wxFile::wxFile::fd_stdout);
  //svgDoc->Save(wxT(outfile));
  svgDoc->Save(wxT("file.svg"));

  // if things bomb, send back a nonzero
  return (nlines);
}

