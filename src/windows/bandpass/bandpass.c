/*
 THIS CODE WAS ADAPTED FOR MAXMSP FROM THE GGEE PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
 */

#include "ext.h"							// standard Max include, always required
#include "ext_obex.h"						// required for new style Max object
#include "m_memory.h"	
#include "filters.h"

////////////////////////// object struct
typedef struct _bandpass
{
	t_object  x_obj;  
    t_float  x_rate;
    t_float  x_freq;
    t_float  x_gain;
    t_float  x_bw;
    
    void *s_float;
} t_bandpass;


///////////////////////// function prototypes
static void bandpass_freq(t_bandpass *x,int f);
static void bandpass_band(t_bandpass *x,int f);
static void bandpass_bang(t_bandpass  *x);
static void bandpass_free(t_bandpass *x);
static void *bandpass_new(t_symbol *s, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *bandpass_class;


int main(void)
{	
	t_class *c;
	
	c = class_new("bandpass", (method)bandpass_new, (method)bandpass_free, (long)sizeof(t_bandpass),0L /* leave NULL!! */,  A_GIMME, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	bandpass_class = c;
    
    addbang((method)bandpass_bang);
    class_addmethod(bandpass_class, (method)bandpass_bang, "bang", 0);
    
    class_addmethod(bandpass_class, (method)bandpass_freq, "int", A_LONG, 0);
    class_addmethod(bandpass_class, (method)bandpass_band, "in1", A_LONG, 0);
    
	post("bandpass version 0.1.0");
	return 0;
}


static void bandpass_bang(t_bandpass *x)
{
    t_atom at[5];
    t_float omega = e_omega(x->x_freq,x->x_rate);
    t_float alpha = e_alpha(x->x_bw* 0.01,omega);
    t_float b1 = 0.;
    t_float b0 = alpha;
    t_float b2 = -alpha;
    t_float a0 = 1 + alpha;
    t_float a1 = -2.*cos(omega);
    t_float a2 = 1 - alpha;
    
        //post("bang %f %f %f",x->x_freq, x->x_gain, x->x_bw); 
    
    if (!check_stability(-a1/a0,-a2/a0,b0/a0,b1/a0,b2/a0)) {
        post("bandpass: filter unstable -> resetting");
        a0=1.;a1=0.;a2=0.;
        b0=1.;b1=0.;b2=0.;
    }
    
    //this is what Max's biquad expects, a bit different from the PD
    // y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2] - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]  
    atom_setfloat(at,b0/a0);
    atom_setfloat(at+1,b1/a0);
    atom_setfloat(at+2,b2/a0);
    atom_setfloat(at+3,a1/a0);
    atom_setfloat(at+4,a2/a0);
    
    outlet_list(x->s_float,0,5,at);

}

static void bandpass_freq(t_bandpass *x,int f)
{
    x->x_freq = f;
    bandpass_bang(x);
}

static void bandpass_band(t_bandpass *x,int f)
{
    x->x_bw = f;
    bandpass_bang(x);
}

static void bandpass_free(t_bandpass *x)
{    
    // free memory    
}

static void *bandpass_new(t_symbol *s, long argc, t_atom *argv)
{
    t_bandpass *x = NULL;
    int i;
    t_atom *ap;

    
    
	if ((x = (t_bandpass *)object_alloc(bandpass_class))) {
        
        // increment ap each time to get to the next atom
        for (i = 0, ap = argv; i < argc; i++, ap++) 
        {
            if (atom_gettype(ap) == A_LONG) 
            {
                if (i==0)  x->x_freq = atom_getlong(ap);
                if (i==1) x->x_bw = atom_getlong(ap);
            }
            else post("%ld: unknown atom type (%ld)", i+1, atom_gettype(ap));
        }        
        
        x->s_float = listout(x);    
        
        intin(x, 1);
        
        x->x_rate = 44100;
    }
    
    return (x);
}

