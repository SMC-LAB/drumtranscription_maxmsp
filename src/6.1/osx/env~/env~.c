/*
  THIS CODE WAS ADAPTED FOR MAXMSP FROM THE TIMBREID PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
 
 bfcc~ - A Bark Frequency Cepstrum Analysis external.
 
 */

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#define t_float float
#define MAXOVERLAP 32
#define INITVSTAKEN 64
#define LOGTEN 2.302585092994

// struct to represent the object's state
typedef struct sigenv {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
	//double			offset; 	// the value of a property of our object
    void *x_outlet;                 /* a "float" outlet */
    void *x_clock;                  /* a "clock" object */
    t_sample *x_buf;                   /* a Hanning window */
    int x_phase;                    /* number of points since last output */
    int x_period;                   /* requested period of output */
    int x_realperiod;               /* period rounded up to vecsize multiple */
    int x_npoints;                  /* analysis window size in samples */
    t_float x_result;                 /* result to output */
    t_sample x_sumbuf[MAXOVERLAP];     /* summing buffer */
    t_float x_f;
    int x_allocforvs;               /* extra buffer for DSP vector size */

} t_sigenv;


// method prototypes
t_float powertodb(t_float f);
void *env_tilde_new(t_symbol *s, long argc, t_atom *argv);
void env_tilde_free(t_sigenv *x);
void env_tilde_dsp(t_sigenv *x, t_signal **sp, short *count);
void env_tilde_dsp64(t_sigenv *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *env_tilde_perform(t_int *w);
void env_tilde_perform64(t_sigenv *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
static void env_tilde_tick(t_sigenv *x);

// global class pointer variable
static t_class *env_tilde_class = NULL;


//***********************************************************************************************

int C74_EXPORT main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("env~", (method)env_tilde_new, (method)env_tilde_free, (long)sizeof(t_sigenv), 0L, A_GIMME, 0);
    
    class_dspinit(c);
	class_register(CLASS_BOX, c);
	env_tilde_class = c;
	
	class_addmethod(c, (method)env_tilde_dsp,		"dsp",		A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 5 and earlier
	class_addmethod(c, (method)env_tilde_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6


	return 0;
}


void *env_tilde_new(t_symbol *s, long argc, t_atom *argv)
{
	t_sigenv *x = (t_sigenv *)object_alloc(env_tilde_class);
    float npoints;
    float period;
    int i, isPow2;
	//s=s;
    
	if (x) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
										// use 0 if you don't need inlets
		//outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
		//x->offset = 0.0;
        
        if(argc > 1)
        {
            npoints = atom_getfloat(argv);  // should perform a check for >64 && power of two
            isPow2 = (int)npoints && !( ((int)npoints-1) & (int)npoints );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead.");
                npoints = 1024;
            }
            else x->x_npoints=(int)npoints;
            
            period = atom_getfloat(argv+1);
            isPow2 = (int)period && !( ((int)period-1) & (int)period );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 512 used instead.");
                period = npoints/2;
            }
            else x->x_period=(int)period;
        }
        else if(argc > 0)
        {
            npoints = atom_getfloat(argv);
            isPow2 = (int)npoints && !( ((int)npoints-1) & (int)npoints );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead.");
                x->x_npoints = 1024;
            }
            else x->x_npoints=(int)npoints;
            
            period = x->x_npoints/2;
        }
        else
        {
            x->x_npoints = 1024;
            x->x_period = 512;
        }
        if (x->x_period < x->x_npoints / MAXOVERLAP + 1)
            x->x_period = x->x_npoints / MAXOVERLAP + 1;
        
      
        if (!(x->x_buf= sysmem_newptrclear(sizeof(t_sample) * (x->x_npoints + INITVSTAKEN))))
        {
            error("env: couldn't allocate buffer");
            return (0);
        }   
        x->x_phase = 0;
        for (i = 0; i < MAXOVERLAP; i++) x->x_sumbuf[i] = 0;
        for (i = 0; i < x->x_npoints; i++)
            x->x_buf[i] = (1. - cos((2 * 3.14159 * i) / x->x_npoints))/x->x_npoints;
        for (; i < x->x_npoints+INITVSTAKEN; i++) x->x_buf[i] = 0;
        
        x->x_clock = clock_new((t_object *)x, (method)env_tilde_tick);
        
        x->x_outlet = floatout(x);
        x->x_f = 0;
        x->x_allocforvs = INITVSTAKEN;
        
              
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void env_tilde_free(t_sigenv *x)
{
    dsp_free((t_pxobject *)x);
    object_free(x->x_clock);
    sysmem_freeptr(x->x_buf);    
}


//***********************************************************************************************
t_float powtodb(t_float f)
{
    if (f <= 0) return (0);
    else
    {
        t_float val = 100 + 10./LOGTEN * log(f);
        return (val < 0 ? 0 : val);
    }
}



//***********************************************************************************************

// this function is called when the DAC is enabled, and "registers" a function for the signal chain in Max 5 and earlier.
// In this case we register the 32-bit, "bfcc_perform" method.
void env_tilde_dsp(t_sigenv *x, t_signal **sp, short *count)
{
	//post("my sample rate is: %f", sp[0]->s_sr);
	
	// dsp_add
	// 1: (t_perfroutine p) perform method
	// 2: (long argc) number of args to your perform method
	// 3...: argc additional arguments, all must be sizeof(pointer) or long
	// these can be whatever, so you might want to include your object pointer in there
	// so that you have access to the info, if you need it.
	dsp_add(env_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
    
    if (x->x_period % sp[0]->s_n) x->x_realperiod =
        x->x_period + sp[0]->s_n - (x->x_period % sp[0]->s_n);
    else x->x_realperiod = x->x_period;
    if (sp[0]->s_n > x->x_allocforvs)
    {
        //void *xx = t_resizebytes_(x->x_buf,
        //                       (x->x_npoints + x->x_allocforvs) * sizeof(t_sample),
        //                       (x->x_npoints + sp[0]->s_n) * sizeof(t_sample));
        if (!(x->x_buf = sysmem_resizeptr(x->x_buf, (x->x_npoints + sp[0]->s_n) * sizeof(t_sample))))
        {
            error("env~: out of memory");
            return;
        }
        //x->x_buf = (t_sample *)xx;
        x->x_allocforvs = sp[0]->s_n;
    }

}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void env_tilde_dsp64(t_sigenv *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	//post("my sample rate is: %f", samplerate);
	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the dsp_add64 arguments are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: a pointer to your object
	// 3: a pointer to your 64-bit perform method
	// 4: flags to alter how the signal chain handles your object -- just pass 0
	// 5: a generic pointer that you can use to pass any additional data to your perform method
	
	object_method(dsp64, gensym("dsp_add64"), x, env_tilde_perform64, 0, NULL);
    
    if (x->x_period % maxvectorsize) x->x_realperiod =
        x->x_period + maxvectorsize - (x->x_period % maxvectorsize);
    else x->x_realperiod = x->x_period;
    
    if (maxvectorsize > x->x_allocforvs)
    {
        //void *xx = t_resizebytes_(x->x_buf,
        //                       (x->x_npoints + x->x_allocforvs) * sizeof(t_sample),
        //                       (x->x_npoints + sp[0]->s_n) * sizeof(t_sample));
        if (!(x->x_buf = sysmem_resizeptr(x->x_buf, (x->x_npoints + maxvectorsize) * sizeof(t_sample))))
        {
            error("env~: out of memory");
            return;
        }
        //x->x_buf = (t_sample *)xx;
        x->x_allocforvs = maxvectorsize;
    }
    

}


//***********************************************************************************************

// this is the 32-bit perform method for Max 5 and earlier
t_int *env_tilde_perform(t_int *w)
{
    
	// DO NOT CALL post IN HERE, but you can call defer_low (not defer)
	
	// args are in a vector, sized as specified in bfcc_dsp method
	// w[0] contains &bfcc_perform, so we start at w[1]
	//t_sigenv *x = (t_sigenv *)(w[1]);
	//t_float *inL = (t_float *)(w[2]);
	//t_float *outL = (t_float *)(w[3]);
	
	t_sigenv *x = (t_sigenv *)(w[1]);
    t_sample *in = (t_float *)(w[2]);
    int n = (int)(w[3]);
    int count;
    t_sample *sump; 
    in += n;
    for (count = x->x_phase, sump = x->x_sumbuf;
         count < x->x_npoints; count += x->x_realperiod, sump++)
    {
        t_sample *hp = x->x_buf + count;
        t_sample *fp = in;
        t_sample sum = *sump;
        int i;
        
        for (i = 0; i < n; i++)
        {
            fp--;
            sum += *hp++ * (*fp * *fp);
        }
        *sump = sum;
    }
    sump[0] = 0;
    x->x_phase -= n;
    if (x->x_phase < 0)
    {
        x->x_result = x->x_sumbuf[0];
        for (count = x->x_realperiod, sump = x->x_sumbuf;
             count < x->x_npoints; count += x->x_realperiod, sump++)
            sump[0] = sump[1];
        sump[0] = 0;
        x->x_phase = x->x_realperiod - n;
        clock_delay(x->x_clock, 0L);
    }
		
	// you have to return the NEXT pointer in the array OR MAX WILL CRASH
	return w + 4;
}


// this is 64-bit perform method for Max 6
void env_tilde_perform64(t_sigenv *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    int n;
	//t_double *inL = ins[0];		// we get audio for each inlet of the object from the **ins argument
    t_sample *in = ins[0];
	//t_double *outL = outs[0];	// we get audio for each outlet of the object from the **outs argument
	n = sampleframes;
    
    
    int count;
    t_sample *sump; 
    in += n;
    for (count = x->x_phase, sump = x->x_sumbuf;
         count < x->x_npoints; count += x->x_realperiod, sump++)
    {
        t_sample *hp = x->x_buf + count;
        t_sample *fp = in;
        t_sample sum = *sump;
        int i;
        
        for (i = 0; i < n; i++)
        {
            fp--;
            sum += *hp++ * (*fp * *fp);
        }
        *sump = sum;
    }
    sump[0] = 0;
    x->x_phase -= n;
    if (x->x_phase < 0)
    {
        x->x_result = x->x_sumbuf[0];
        for (count = x->x_realperiod, sump = x->x_sumbuf;
             count < x->x_npoints; count += x->x_realperiod, sump++)
            sump[0] = sump[1];
        sump[0] = 0;
        x->x_phase = x->x_realperiod - n;
        clock_delay(x->x_clock, 0L);
    }
    
	// this perform method simply copies the input to the output, offsetting the value
	//while (n--)
	//	*outL++ = *inL++ + x->offset;
}

static void env_tilde_tick(t_sigenv *x) /* callback function for the clock */
{
    outlet_float(x->x_outlet, powtodb(x->x_result));
}

