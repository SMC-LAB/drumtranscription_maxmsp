/* Copyright (c) 1997-1999 Miller Puckette.
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */
#include "m_memory.h"


void *t_getbytes_(long nbytes)
{
	return sysmem_newptr(nbytes);
}

void t_freebytes_(void *fatso, long nbytes)
{
	sysmem_freeptr (fatso);
}
void *t_resizebytes_(void *old, long oldsize, long newsize)
{
	return sysmem_resizeptr(old, newsize);
}