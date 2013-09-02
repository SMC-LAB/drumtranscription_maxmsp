/* Copyright (c) 1997-1999 Miller Puckette.
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */
#include <stdlib.h>
#include <string.h>
#include "ext.h"	

/* --------------- memory management -------------------- */
void *t_getbytes_(long nbytes);
//void *t_getzbytes_(size_t nbytes);
//void *t_copybytes_(void *src, long nbytes);
void t_freebytes_(void *fatso, long nbytes);
void *t_resizebytes_(void *old, long oldsize, long newsize);