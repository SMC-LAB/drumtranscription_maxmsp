/**
   THIS CODE WAS ADAPTED FOR MAXMSP FROM THE TIMBREID PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
	@file
	specFlatness~: a simple audio object for Max
	original by: jeremy bernstein, jeremy@bootsquad.com	
	@ingroup examples	
*/

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include "tIDLib.h"
#define t_float float

// struct to represent the object's state
typedef struct _specFlatness {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
	//double			offset; 	// the value of a property of our object
    t_float sr;
    t_float n;
	int windowFunction;
    int overlap;
    int powerSpectrum;
    int window;
    double lastDspTime;
    t_sample *signal_R;
    t_float *nthRoots;
    t_float *blackman;
    t_float *cosine;
    t_float *hamming;
    t_float *hann;    
    t_float x_f;
    
    void *x_flatness;

} t_specFlatness;


// method prototypes
static void specFlatness_tilde_bang(t_specFlatness *x);
static void specFlatness_tilde_window(t_specFlatness *x, int w);
static void specFlatness_tilde_overlap(t_specFlatness *x, int o);
static void specFlatness_tilde_windowFunction(t_specFlatness *x, int f);
static void specFlatness_tilde_powerSpectrum(t_specFlatness *x, int spec);

void *specFlatness_new(t_symbol *s, long argc, t_atom *argv);
void specFlatness_free(t_specFlatness *x);
void specFlatness_dsp(t_specFlatness *x, t_signal **sp, short *count);
void specFlatness_dsp64(t_specFlatness *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *specFlatness_perform(t_int *w);
void specFlatness_perform64(t_specFlatness *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *specFlatness_class = NULL;


//***********************************************************************************************

int main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("specFlatness~", (method)specFlatness_new, (method)dsp_free, (long)sizeof(t_specFlatness), 0L, A_GIMME, 0);
    
    class_dspinit(c);
	class_register(CLASS_BOX, c);
	specFlatness_class = c;
	
	class_addmethod(c, (method)specFlatness_dsp,		"dsp",		A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 5 and earlier
	class_addmethod(c, (method)specFlatness_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
	//class_addmethod(c, (method)specFlatness_assist,	"assist",	A_CANT, 0);
    
    
	addbang((method)specFlatness_tilde_bang);
    
    class_addmethod(specFlatness_class, 
                    (method)specFlatness_tilde_bang,
                    "bang",
                    0
                    ); 
    
    
	class_addmethod(
                    specFlatness_class,
                    (method)specFlatness_tilde_window,
                    "window",
                    A_LONG,
                    0
                    );
    
	class_addmethod(
                    specFlatness_class,
                    (method)specFlatness_tilde_overlap,
                    "overlap",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    specFlatness_class,
                    (method)specFlatness_tilde_windowFunction,
                    "window_function",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    specFlatness_class,
                    (method)specFlatness_tilde_powerSpectrum,
                    "power_spectrum",
                    A_DEFLONG,
                    0
                    );
	


	return 0;
}


void *specFlatness_new(t_symbol *s, long argc, t_atom *argv)
{
	t_specFlatness *x = (t_specFlatness *)object_alloc(specFlatness_class);
    int i, isPow2;
    double timef = 0.0;
	//s=s;
    
	if (x) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
										// use 0 if you don't need inlets
		//outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
		//x->offset = 0.0;
        
        
        x->x_flatness = floatout(x);
        //x->x_flatness = outlet_new(&x->x_obj, &s_float);        
        
        if(argc > 1)
        {
            x->window = atom_getfloat(argv);
            isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead");
                x->window = 1024;
            };
                      
        }
        else if(argc > 0)
        {
            x->window = atom_getfloat(argv);
            isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead");
                x->window = 1024;
            };
            
            
        }
        else
        {
            x->window = 1024;
        };
        
        x->sr = 44100.0;
        x->n = 64.0;
        x->overlap = 1;
        x->windowFunction = 4; // 4 is hann window
        x->powerSpectrum = 0;
        //x->lastDspTime = clock_getlogicaltime();
        clock_getftime(&timef);
        x->lastDspTime = timef;
        
        x->signal_R = (t_sample *)t_getbytes_((x->window+x->n) * sizeof(t_sample));
        x->nthRoots = (t_float *)t_getbytes((x->window*0.5+1)*sizeof(t_float));
        
        for(i=0; i<(x->window+x->n); i++)
            x->signal_R[i] = 0.0;
        
        x->blackman = (t_float *)t_getbytes_(x->window*sizeof(t_float));
        x->cosine = (t_float *)t_getbytes_(x->window*sizeof(t_float));
        x->hamming = (t_float *)t_getbytes_(x->window*sizeof(t_float));
        x->hann = (t_float *)t_getbytes_(x->window*sizeof(t_float));
        
        // initialize signal windowing functions
        tIDLib_blackmanWindow(x->blackman, x->window);
        tIDLib_cosineWindow(x->cosine, x->window);
        tIDLib_hammingWindow(x->hamming, x->window);
        tIDLib_hannWindow(x->hann, x->window);
        
        post("specFlatness~: window size: %i", (int)x->window);
        
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void specFlatness_free(t_specFlatness *x) 
{
	// free the input buffer memory
    t_freebytes_(x->signal_R, (x->window+x->n)*sizeof(t_sample));
    
	t_freebytes(x->nthRoots, (x->window*0.5+1)*sizeof(t_float));
    
	// free the window memory
	t_freebytes_(x->blackman, x->window*sizeof(t_float));
	t_freebytes_(x->cosine, x->window*sizeof(t_float));
	t_freebytes_(x->hamming, x->window*sizeof(t_float));
	t_freebytes_(x->hann, x->window*sizeof(t_float));
}


//***********************************************************************************************

static void specFlatness_tilde_bang(t_specFlatness *x)
{
    int i, j, window, windowHalf, bangSample;
    t_sample *signal_R, *signal_I;
    t_float dividend, divisor, windowHalfPlusOneRecip, flatness, *windowFuncPtr;
	double currentTime,timef =0.0;
    
    window = x->window;
    windowHalf = window*0.5;
    windowHalfPlusOneRecip = 1.0/(t_float)(windowHalf+1);
	
	// create local memory
   	signal_R = (t_sample *)t_getbytes(window*sizeof(t_sample));
	signal_I = (t_sample *)t_getbytes((windowHalf+1)*sizeof(t_sample));
    
	//currentTime = clock_gettimesince(x->lastDspTime);
    clock_getftime(&timef);
    currentTime =  timef - x->lastDspTime;
    
	bangSample = (int)(((currentTime/1000.0)*x->sr)+0.5); // round
    
	if (bangSample < 0)
        bangSample = 0;
    else if ( bangSample >= x->n )
        bangSample = x->n - 1;
    
	// construct analysis window using bangSample as the end of the window
	for(i=0, j=bangSample; i<window; i++, j++)
		signal_R[i] = x->signal_R[j];
    
	// set window function
	windowFuncPtr = x->hann; //default case to get rid of compile warning
    
	switch(x->windowFunction)
	{
		case 0:
			break;
		case 1:
			windowFuncPtr = x->blackman;
			break;
		case 2:
			windowFuncPtr = x->cosine;
			break;
		case 3:
			windowFuncPtr = x->hamming;
			break;
		case 4:
			windowFuncPtr = x->hann;
			break;
		default:
			break;
	};
    
	// if windowFunction == 0, skip the windowing (rectangular)
	if(x->windowFunction>0)
		for(i=0; i<window; i++, windowFuncPtr++)
			signal_R[i] *= *windowFuncPtr;
    
	mayer_realfft(window, signal_R);
	tIDLib_realfftUnpack(window, windowHalf, signal_R, signal_I);
	tIDLib_power(windowHalf+1, signal_R, signal_I);
    
	if(!x->powerSpectrum)
		tIDLib_mag(windowHalf+1, signal_R);
    
	dividend=1; // to get the product of all terms for geometric mean
	divisor=0;
	flatness=0;
	
	// geometric mean
	// take the nth roots first so as not to lose data.
	for(i=0; i<=windowHalf; i++)
        x->nthRoots[i] = pow(signal_R[i], windowHalfPlusOneRecip);
    
	// take the product of nth roots
	// what to do with values that are zero?  for now, ignoring them.	
	for(i=0; i<=windowHalf; i++)
		if(x->nthRoots[i] != 0)
			dividend *= x->nthRoots[i];
    
	for(i=0; i<=windowHalf; i++)
		divisor += signal_R[i];
	
	divisor *= windowHalfPlusOneRecip; // arithmetic mean
	
	if(divisor==0)
		divisor = 1;
	
	flatness = dividend/divisor;
	
	outlet_float(x->x_flatness, flatness);
    
	// free local memory
	t_freebytes(signal_R, window*sizeof(t_sample));
	t_freebytes(signal_I, (windowHalf+1)*sizeof(t_sample));
}


static void specFlatness_tilde_window(t_specFlatness *x, int w)
{
    int i, window, isPow2;
	
	window = w;
	
	isPow2 = window && !( (window-1) & window );
	
	if( !isPow2 )
		error("requested window size is not a power of 2");
	else
	{
		x->signal_R = (t_sample *)t_resizebytes_(x->signal_R, (x->window+x->n) * sizeof(t_sample), (window+x->n) * sizeof(t_sample));
		x->blackman = (t_float *)t_resizebytes_(x->blackman, x->window*sizeof(t_float), window*sizeof(t_float));
		x->cosine = (t_float *)t_resizebytes_(x->cosine, x->window*sizeof(t_float), window*sizeof(t_float));
		x->hamming = (t_float *)t_resizebytes_(x->hamming, x->window*sizeof(t_float), window*sizeof(t_float));
		x->hann = (t_float *)t_resizebytes_(x->hann, x->window*sizeof(t_float), window*sizeof(t_float));
		x->nthRoots = (t_float *)t_resizebytes_(x->nthRoots, (x->window*0.5+1)*sizeof(t_float), (window*0.5+1)*sizeof(t_float));
        
		x->window = (t_float)window;
        
		// re-init window functions
		tIDLib_blackmanWindow(x->blackman, x->window);
		tIDLib_cosineWindow(x->cosine, x->window);
		tIDLib_hammingWindow(x->hamming, x->window);
		tIDLib_hannWindow(x->hann, x->window);
        
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
        
		for(i=0; i<=x->window*0.5; i++)
			x->nthRoots[i] = 0.0;
		
		post("window size: %i", (int)x->window);
	}

}


static void specFlatness_tilde_overlap(t_specFlatness *x, int o)
{
	int overlap;
    
	overlap = o;
    
	// this change will be picked up in _dsp, where things will be updated based on the samplerate sp[0]->s_sr/x->overlap;
	if(overlap > 0)
		x->overlap = overlap;
	else
		error("overlap must be at least 1.");
    
    post("overlap: %i", x->overlap);
}


static void specFlatness_tilde_windowFunction(t_specFlatness *x, int f)
{
    f = (f<0)?0:f;
    f = (f>4)?4:f;
	x->windowFunction = f;
    
	switch(x->windowFunction)
	{
		case 0:
			post("window function: rectangular.");
			break;
		case 1:
			post("window function: blackman.");
			break;
		case 2:
			post("window function: cosine.");
			break;
		case 3:
			post("window function: hamming.");
			break;
		case 4:
			post("window function: hann.");
			break;
		default:
			break;
	};
}


// magnitude spectrum == 0, power spectrum == 1
static void specFlatness_tilde_powerSpectrum(t_specFlatness *x, int spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->powerSpectrum = spec;
    
	if(x->powerSpectrum)
		post("using power spectrum");
	else
		post("using magnitude spectrum");
}



//***********************************************************************************************

// this function is called when the DAC is enabled, and "registers" a function for the signal chain in Max 5 and earlier.
// In this case we register the 32-bit, "specFlatness_perform" method.
void specFlatness_dsp(t_specFlatness *x, t_signal **sp, short *count)
{
    int i;
    double timef;
	post("my sample rate is: %f", sp[0]->s_sr);
	
	// dsp_add
	// 1: (t_perfroutine p) perform method
	// 2: (long argc) number of args to your perform method
	// 3...: argc additional arguments, all must be sizeof(pointer) or long
	// these can be whatever, so you might want to include your object pointer in there
	// so that you have access to the info, if you need it.
	dsp_add(specFlatness_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
    
    
    // compare n to stored n and recalculate filterbank if different
	if(sp[0]->s_sr != (x->sr*x->overlap) || sp[0]->s_n != x->n)
	{
		x->signal_R = (t_sample *)t_resizebytes_(x->signal_R, (x->window+x->n) * sizeof(t_sample), (x->window+sp[0]->s_n) * sizeof(t_sample));
        
		x->sr = sp[0]->s_sr/x->overlap;
		x->n = sp[0]->s_n;
		//x->lastDspTime = clock_getlogicaltime();
        clock_getftime(&timef);
        x->lastDspTime =  timef;
        
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
                
    	post("specFlatness~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)x->sr, (int)x->n);
	};
}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void specFlatness_dsp64(t_specFlatness *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    int i;
    double timef;
	post("my sample rate is: %f", samplerate);
	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the dsp_add64 arguments are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: a pointer to your object
	// 3: a pointer to your 64-bit perform method
	// 4: flags to alter how the signal chain handles your object -- just pass 0
	// 5: a generic pointer that you can use to pass any additional data to your perform method
	
	object_method(dsp64, gensym("dsp_add64"), x, specFlatness_perform64, 0, NULL);
    
    
    // compare n to stored n and recalculate filterbank if different
	if(samplerate != (x->sr*x->overlap) || maxvectorsize != x->n)
	{
		x->signal_R = (t_sample *)t_resizebytes_(x->signal_R, (x->window+x->n) * sizeof(t_sample), (x->window+maxvectorsize) * sizeof(t_sample));
        
		x->sr = samplerate/x->overlap;
		x->n = maxvectorsize;
        
		//x->lastDspTime = clock_getlogicaltime();
        clock_getftime(&timef);
        x->lastDspTime =  timef;
        
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
        
    	post("specFlatness~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)x->sr, (int)x->n);
	};
}


//***********************************************************************************************

// this is the 32-bit perform method for Max 5 and earlier
t_int *specFlatness_perform(t_int *w)
{
    int i, n;
    double timef;
	// DO NOT CALL post IN HERE, but you can call defer_low (not defer)
	
	// args are in a vector, sized as specified in specFlatness_dsp method
	// w[0] contains &specFlatness_perform, so we start at w[1]
	t_specFlatness *x = (t_specFlatness *)(w[1]);
	//t_float *inL = (t_float *)(w[2]);
	//t_float *outL = (t_float *)(w[3]);
    t_sample *in = (t_float *)(w[2]);
	n = (int)w[3];
	
	// shift signal buffer contents back.
	for(i=0; i<x->window; i++)
		x->signal_R[i] = x->signal_R[i+n];
    
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)x->window+i] = in[i];
    
	//x->lastDspTime = clock_getlogicaltime();
    clock_getftime(&timef);
    x->lastDspTime =  timef;
		
	// you have to return the NEXT pointer in the array OR MAX WILL CRASH
	return w + 4;
}


// this is 64-bit perform method for Max 6
void specFlatness_perform64(t_specFlatness *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    int i, n;
    double timef;
	//t_double *inL = ins[0];		// we get audio for each inlet of the object from the **ins argument
    t_sample *in = ins[0];
	//t_double *outL = outs[0];	// we get audio for each outlet of the object from the **outs argument
	n = sampleframes;
    
    // shift signal buffer contents back.
	for(i=0; i<x->window; i++)
		x->signal_R[i] = x->signal_R[i+n];
    
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)x->window+i] = in[i];
    
	//x->lastDspTime = clock_getlogicaltime();
    clock_getftime(&timef);
    x->lastDspTime =  timef;

	
	// this perform method simply copies the input to the output, offsetting the value
	//while (n--)
	//	*outL++ = *inL++ + x->offset;
}

