
/* onplot_movie */

#include "gprim.h"

int colormap;
int backgnd;			// background color
int textcol;			// text color
int erasecol;			// erase color
int make_movie;
int show_inputs;
int space_time;
int volt_map;
int inact_map;
int volt_inact;
int dend_layers;
int inputs_backgr;		// draw inputs before dendrites
int noplots;
int dend_vcolor1;
int dend_vcolor2;
int denddens_color1;
int denddens_color1t;
int denddens_color2;
int denddens_color2t;
const char *movie_title1 = NULL;
const char *movie_title2 = NULL;

extern int rec_ct;
extern int rec_cn;

extern double Vmax;
extern double Vmin;

/*------------------------------------------------*/

void defparams_onplot_movie(void)

{
  setptr("make_movie",     &make_movie);
  setptr("show_inputs",    &show_inputs);
  setptr("space_time",     &space_time);
  setptr("volt_map",       &volt_map);
  setptr("inact_map",      &inact_map);
  setptr("volt_inact",     &volt_inact);
  setptr("dend_layers",    &dend_layers);
  setptr("inputs_backgr",  &inputs_backgr);
  setptr("colormap",	   &colormap);
  setptr("backgnd",	   &backgnd);
  setptr("noplots",	   &noplots);
  setptr("dend_vcolor1",   &dend_vcolor1);
  setptr("dend_vcolor2",   &dend_vcolor2);
  setptr("denddens_color1", &denddens_color1);
  setptr("denddens_color1t",&denddens_color1t);
  setptr("denddens_color2", &denddens_color2);
  setptr("denddens_color2t", &denddens_color2t);
}

/*------------------------------------------------*/

static double oldtime;

void draw_time(double dtime)

{
   char sbuf[80];

  gframe ("time");
  if (dtime < 1e-6 && dtime > -1e-6) dtime = 0;
  gcwidth(0.025);
  gmove(0,0);
  gpen (erasecol);		// erase letters
  sprintf(sbuf,"%-7.4f s",oldtime);
  gtext (sbuf);
  gmove(0,0);
  gpen (textcol);		// write white on black, black on white
  sprintf(sbuf,"%-7.4f s",dtime);
  gtext (sbuf);
  oldtime = dtime;
  gframe ("..");
}

/*------------------------------------------------*/

void drawvcolscale(const char *frame, double cmin, double cmax, const char *unit, int colormap)

 /* Draw the scalebar for cacolor display */

{
    int i,colbase,numcols;
    double colbarlen,dist,wid,labellen;
    double x1,y1,x2,y2,x3,y3,x4,y4;
    char nbuf[80];
    static int colors[10][2] = {0,16,16,16,32,32,64,16,80,16,96,16,112,16,128,100,228,100,0,0}; 
				/* see colormap in manual */
  gframe (frame);
  gcwidth(0.02);
  colbarlen=0.36;
  wid = 0.03;
  colbase = colors[colormap][0];	/* lowest color in colormap */
  numcols = colors[colormap][1];	/* number of colors used */
  dist = colbarlen/numcols;

  for (i=0; i<numcols; i++){
        gpen(i+colbase); 
	x1 = wid/2;
	y1 = i*dist;
	x2 = x1;
	y2 = (i+1)*dist;
	x3 = -wid/2;
	y3 = y2;
	x4 = x3;
	y4 = y1;
	grect(x1,y1,x2,y2,x3,y3,x4,y4,1);
  }

  gpen (textcol);
  sprintf (nbuf,"%-.3g %s",cmax,unit);        /* highest value */
  labellen = strlen(nbuf);
  gmove(-labellen/2*0.015,colbarlen+0.01);
  gtext(nbuf);

  sprintf (nbuf,"%-.3g %s",cmin,unit);          /* lowest value  */
  labellen = strlen(nbuf);
  gmove(-labellen/2*0.015,-0.03);
  gtext(nbuf);
  gframe ("..");
}

/*------------------------------------------------*/

void onplot_movie_init(void)

{
  if (notinit(make_movie))   make_movie = 0;
  if (notinit(show_inputs)) show_inputs = 0;
  if (notinit(space_time))   space_time = 0;
  if (notinit(volt_map))       volt_map = 1;
  if (notinit(inact_map))     inact_map = 0;
  if (volt_map && inact_map) volt_inact = 1;
  else                       volt_inact = 0;
  if (volt_inact)            space_time = 0;
  if (notinit(dend_layers)) dend_layers = 0;
  // if (dend_layers)	       volt_map = 0;

  if (notinit(inputs_backgr))     inputs_backgr = 1;	// draw inputs before dendrites
  if (notinit(denddens_color1)) denddens_color1 = _COL;	// color row in density file
  if (notinit(denddens_color2)) denddens_color2 = _COL;
  if (notinit(denddens_color1t)) denddens_color1t = _COL;  // color row in density file for labels
  if (notinit(denddens_color2t)) denddens_color2t = _COL;

  if (notinit(backgnd))         backgnd = 0;		// background color, need for text
  if (backgnd==0) textcol = 7;				// secure black (0 doesn't print)
  else            textcol = GCOL1;			// use secure black (0 doesn't print) 
  erasecol = (textcol==7 ? GCOL1 : GCOL99);		// use secure white (white swaps to black in vid.c)

  if (make_movie && !space_time && !dend_layers) noplots = 1;
  if (notinit(noplots)) noplots = 0;

  if (make_movie) {	 /* set up movie and picture components */

   disp_calib_len = 0;      // turn off regular calib line at bottom
   oldtime = 0;
   if (notinit(colormap)) colormap = 2;	/* set colormap for dsgc */
   if (space_time) {
     gframe ("/movie");
     gorigin (0.15,0.37);	/* movie smaller, slightly up */
     gsize (0.68);
     gmove (0.3,0.9);
     gpen (textcol);		// write white on black, black on white
     gtext (movie_title1);
     gframe ("time");
     gorigin(0.50,0.015);	/* set time position */
     gframe ("..");
     gframe ("..");
     gframe ("/gc_col_bar");
     gorigin (0.90,0.45);	/* set bar position */
     gframe ("..");
   } else if (volt_inact || dend_layers) {/* dual disp volts, na inact; or 2 dend layers */
     gframe ("/movie");
     gorigin (0.15,0.55);		/* movie smaller, halfway up */
     gsize (0.5);
     gmove (0.05,0.85);
     gpen (textcol);		// write white on black, black on white
     gtext (movie_title1);
     gframe ("time");
     gorigin(0.60,0.85);	/* set time position */
     gframe ("..");
     gframe ("..");
     gframe ("/movie2");
     gorigin (0.15,0.25);	/* movie smaller, bottom */
     gsize (0.5);
     gmove (0.05,0.8);
     gpen (textcol);		// write white on black, black on white
     gtext (movie_title2);
     gframe ("..");
     gframe ("/gc_col_bar");
     gorigin (0.90,0.40);	/* set bar position */
     gframe ("..");
   } else {
     gframe ("/movie");
     gorigin (0,0);		/* movie fills screen */
     gmove (0.3,0.95);
     gpen (textcol);		// write white on black, black on white
     gtext (movie_title1);
     gframe ("time");
     gorigin(0.55,0.015);	/* set time position */
     gframe ("..");
     gframe ("..");
     gframe ("/gc_col_bar");
     gorigin (0.94,0.20);	/* set bar position */
     gframe ("..");
   } 
   if       (volt_map) drawvcolscale("gc_col_bar",Vmin,Vmax,"V",colormap);
   else if (inact_map) drawvcolscale("gc_col_bar",0,0.5,"F",colormap);
  fflush(stdout);
  }
}

/*------------------------------------------------*/

void onplot_movie(void) 

{
   int cmap, only, color, voltmap_color1, voltmap_color2;
   double vmax, vmin, dscale;
   static int runyet = 0;

  if (make_movie) {
    if (!notinit(dend_vcolor1)) voltmap_color1  = dend_vcolor1;
    else		        voltmap_color1  = vcolor;
    if (!notinit(dend_vcolor2)) voltmap_color2 = dend_vcolor2;
    else		        voltmap_color2 = vcolor;
    if (simtime >= setxmin && (int(simtime/ploti+0.0001) % int((frame_int/ploti))) == 0) {
      disp |= DMOVIE;		/* turn on movie display */
      if (!runyet) {
         gframe ("/movie");
	 gcwidth (0.025);
	 if (dend_layers)
	   disp_calib(0.95,0.23,set_int_val(dispsize/5),dispsize,cyan);
	 else
	   disp_calib(0.95,0.08,set_int_val(dispsize/5),dispsize,cyan);
	 runyet = 1;
         gframe ("..");
      }
      if (volt_inact) {
        gframe ("/movie");
        draw_time(simtime);
        if (inputs_backgr && show_inputs) draw_inputs();
        set_rcolors(rec_ct,rec_cn);                /* set region colors */
        display(MATCHING,ndt(rec_ct,rec_cn,-1), only=1, color=voltmap_color1,
				  vmax=Vmax,vmin=Vmin,cmap=colormap,dscale=1);
        if (!inputs_backgr && show_inputs) draw_inputs();
        gframe ("..");
        gframe ("/movie2");		/* show inactivation */
        display (MATCHING,ndt(rec_ct,rec_cn,-1), only=1,color=naicolor,
				  vmax=0.5,vmin=0, cmap=colormap,dscale=1);
        gframe ("..");
      }

      else if (dend_layers) {	// make movie of 2 separate dendritic layers
        gframe ("/movie");		// show first layer
        draw_time(simtime);
	_COL = denddens_color1;		// set colors in density file 
        set_rcolors(rec_ct,rec_cn);                /* set region colors */
        if (inputs_backgr && show_inputs) draw_inputs();
        display(MATCHING,ndt(rec_ct,rec_cn,-1), only=1, color=voltmap_color1,
				  vmax=Vmax,vmin=Vmin,cmap=colormap,dscale=1);
	_COL = denddens_color1t;		// set colors in density file 
        set_rcolors(rec_ct,rec_cn);                /* set region colors */
         display (LABEL, MATCHING,  ndt(rec_ct,rec_cn,-1), node_color=RCOLOR, dsgc_nscale);
        if (!inputs_backgr && show_inputs) draw_inputs();
        gframe ("..");
        gframe ("/movie2");		/* show second layer */
	_COL = denddens_color2;		// set colors in density file 
        set_rcolors(rec_ct,rec_cn);                /* set region colors */
        display (MATCHING,ndt(rec_ct,rec_cn,-1), only=1,color=voltmap_color2,
				  vmax=Vmax,vmin=Vmin, cmap=colormap,dscale=1);
	_COL = denddens_color2t;		// set colors in density file 
        set_rcolors(rec_ct,rec_cn);                /* set region colors */
        display (LABEL, MATCHING,  ndt(rec_ct,rec_cn,-1), node_color=RCOLOR, dsgc_nscale);
        gframe ("..");
      }

      else {
       gframe ("/movie");
       draw_time(simtime);
       set_rcolors(rec_ct,rec_cn);                /* set region colors */
       if (inputs_backgr && show_inputs) draw_inputs();
       if (volt_map) display (MATCHING,ndt(rec_ct,rec_cn,-1), only=1,color=voltmap_color1,
				  vmax=Vmax,vmin=Vmin,cmap=colormap,dscale=1);
       if (inact_map) display(MATCHING,ndt(rec_ct,rec_cn,-1), only=1,color=naicolor,
				  vmax=0.5,vmin=0,cmap=colormap,dscale=1);
        display (LABEL, MATCHING,  ndt(rec_ct,rec_cn,-1), node_color=RCOLOR, dsgc_nscale);
       if (!inputs_backgr && show_inputs) draw_inputs();

      //display stim dscale 5;      /* show light stimulus */
      gframe ("..");
     };
     display_page();
     disp &= ~DMOVIE;		/* turn off movie display */
    }
  }
}
