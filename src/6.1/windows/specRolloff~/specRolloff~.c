/**
  THIS CODE WAS ADAPTED FOR MAXMSP FROM THE TIMBREID PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
	@file
	specRolloff~: a simple audio object for Max
	original by: jeremy bernstein, jeremy@bootsquad.com	
	@ingroup examples	
*/

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include "tIDLib.h"
#define t_float float

// struct to represent the object's state
typedef struct _specRolloff {
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
    t_float *binFreqs;
    t_float concentration;
    t_float *blackman;
    t_float *cosine;
    t_float *hamming;
    t_float *hann;    
    t_float x_f;
    
    void *x_rolloff;

} t_specRolloff;


// method prototypes
static void specRolloff_tilde_bang(t_specRolloff *x);
static void specRolloff_tilde_concentration(t_specRolloff *x, t_float c);
static void specRolloff_tilde_window(t_specRolloff *x, int w);
static void specRolloff_tilde_overlap(t_specRolloff *x, int o);
static void specRolloff_tilde_windowFunction(t_specRolloff *x, int f);
static void specRolloff_tilde_powerSpectrum(t_specRolloff *x, int spec);

void *specRolloff_new(t_symbol *s, long argc, t_atom *argv);
void specRolloff_free(t_specRolloff *x);
void specRolloff_dsp(t_specRolloff *x, t_signal **sp, short *count);
void specRolloff_dsp64(t_specRolloff *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *specRolloff_perform(t_int *w);
void specRolloff_perform64(t_specRolloff *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *specRolloff_class = NULL;


//***********************************************************************************************

int C74_EXPORT main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("specRolloff~", (method)specRolloff_new, (method)dsp_free, (long)sizeof(t_specRolloff), 0L, A_GIMME, 0);
    
    class_dspinit(c);
	class_register(CLASS_BOX, c);
	specRolloff_class = c;
	
	class_addmethod(c, (method)specRolloff_dsp,		"dsp",		A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 5 and earlier
	class_addmethod(c, (method)specRolloff_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
	//class_addmethod(c, (method)specRolloff_assist,	"assist",	A_CANT, 0);
    
    
	addbang((method)specRolloff_tilde_bang);
    class_addmethod(specRolloff_class, 
                    (method)specRolloff_tilde_bang,
                    "bang",
                    0
                    ); 
    
	class_addmethod(
                    specRolloff_class,
                    (method)specRolloff_tilde_concentration,
                    "concentration",
                    A_DEFFLOAT,
                    0
                    );
    
	class_addmethod(
                    specRolloff_class,
                    (method)specRolloff_tilde_window,
                    "window",
                    A_LONG,
                    0
                    );
    
	class_addmethod(
                    specRolloff_class,
                    (method)specRolloff_tilde_overlap,
                    "overlap",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    specRolloff_class,
                    (method)specRolloff_tilde_windowFunction,
                    "window_function",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    specRolloff_class,
                    (method)specRolloff_tilde_powerSpectrum,
                    "power_spectrum",
                    A_DEFLONG,
                    0
                    );
	


	return 0;
}


void *specRolloff_new(t_symbol *s, long argc, t_atom *argv)
{
	t_specRolloff *x = (t_specRolloff *)object_alloc(specRolloff_class);
    int i, isPow2;
    double timef = 0.0;
	//s=s;
    
	if (x) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
										// use 0 if you don't need inlets
		//outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
		//x->offset = 0.0;
        
        
        x->x_rolloff = floatout(x);
        //x->x_rolloff = outlet_new(&x->x_obj, &s_float);        
        
        if(argc > 1)
        {
            x->window = atom_getfloat(argv);
            isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead");
                x->window = 1024;
            };
            
            x->concentration = atom_getfloat(argv+1);
            if(x->concentration<0)
            {
                error("concentration must be between 0.0 and 1.0.");
                x->concentration = 0.85;
            }
            else if(x->concentration>1.0)
            {
                error("concentration must be between 0.0 and 1.0.");
                x->concentration = 0.85;
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
            
            x->concentration = 0.85;
            
        }
        else
        {
            x->window = 1024;
            x->concentration = 0.85;
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
        
        x->binFreqs = (t_float *)t_getbytes_(x->window*sizeof(t_float));
        
        // freqs for each bin based on current window size and sample rate
        for(i=0; i<x->window; i++)
            x->binFreqs[i] = (x->sr/x->window)*i;
        
        post("specRolloff~: window size: %i", (int)x->window);
        
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void specRolloff_free(t_specRolloff *x) 
{
	// free the input buffer memory
    t_freebytes_(x->signal_R, (x->window+x->n)*sizeof(t_sample));
    
	// free the binFreqs memory
    t_freebytes_(x->binFreqs, x->window*sizeof(t_sample));
    
	// free the window memory
	t_freebytes_(x->blackman, x->window*sizeof(t_float));
	t_freebytes_(x->cosine, x->window*sizeof(t_float));
	t_freebytes_(x->hamming, x->window*sizeof(t_float));
	t_freebytes_(x->hann, x->window*sizeof(t_float));
}


//***********************************************************************************************

static void specRolloff_tilde_bang(t_specRolloff *x)
{
    int i, j, window, windowHalf, bangSample;
    t_sample *signal_R, *signal_I;
    t_float energy, energyTarget, rolloff, *windowFuncPtr;
	double currentTime,timef;
    
    window = x->window;
    windowHalf = window*0.5;
    
	// create local memory
   	signal_R = (t_sample *)t_getbytes_(window*sizeof(t_sample));
	signal_I = (t_sample *)t_getbytes_((windowHalf+1)*sizeof(t_sample));
    
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
    
	energyTarget=0;
    
	for(i=0; i<=windowHalf; i++)
		energyTarget += signal_R[i];
    
	energyTarget *= x->concentration;
    
	energy=0;
	i=0;
    
	while(energy <= energyTarget)
	{
		energy += signal_R[i];
		i++;
        
		if(i>windowHalf)
			break;
	}
    
	i = (i==0)?1:i;
    
	rolloff = x->binFreqs[i-1]; // back up one because the last one went over...
    
	outlet_float(x->x_rolloff, rolloff);
    
	// free local memory
	t_freebytes_(signal_R, window*sizeof(t_sample));
	t_freebytes_(signal_I, (windowHalf+1)*sizeof(t_sample));
}


static void specRolloff_tilde_concentration(t_specRolloff *x, t_float c)
{
	if(c<0 || c>1.0)
		error("concentration must be between 0.0 and 1.0.");
	else
		x->concentration = c;
    
    post("concentration: %f", x->concentration);
    
}


static void specRolloff_tilde_window(t_specRolloff *x, int w)
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
		x->binFreqs = (t_float *)t_resizebytes_(x->binFreqs, x->window*sizeof(t_float), window*sizeof(t_float));
        
		x->window = (t_float)window;
        
		// re-init window functions
		tIDLib_blackmanWindow(x->blackman, x->window);
		tIDLib_cosineWindow(x->cosine, x->window);
		tIDLib_hammingWindow(x->hamming, x->window);
		tIDLib_hannWindow(x->hann, x->window);
        
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
        
		// freqs for each bin based on current window size and sample rate
		for(i=0; i<x->window; i++)
			x->binFreqs[i] = (x->sr/x->window)*i;
        
		post("window size: %i", (int)x->window);
	}
}


static void specRolloff_tilde_overlap(t_specRolloff *x, int o)
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


static void specRolloff_tilde_windowFunction(t_specRolloff *x, int f)
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
static void specRolloff_tilde_powerSpectrum(t_specRolloff *x, int spec)
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
// In this case we register the 32-bit, "specRolloff_perform" method.
void specRolloff_dsp(t_specRolloff *x, t_signal **sp, short *count)
{
    int i;
    double timef;
	//post("my sample rate is: %f", sp[0]->s_sr);
	
	// dsp_add
	// 1: (t_perfroutine p) perform method
	// 2: (long argc) number of args to your perform method
	// 3...: argc additional arguments, all must be sizeof(pointer) or long
	// these can be whatever, so you might want to include your object pointer in there
	// so that you have access to the info, if you need it.
	dsp_add(specRolloff_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
    
    
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
        
		// freqs for each bin based on current window size and sample rate
		for(i=0; i<x->window; i++)
			x->binFreqs[i] = (x->sr/x->window)*i;
        
    	post("specRolloff~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)x->sr, (int)x->n);
	};
}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void specRolloff_dsp64(t_specRolloff *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    int i;
    double timef;
	//post("my sample rate is: %f", samplerate);
	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the dsp_add64 arguments are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: a pointer to your object
	// 3: a pointer to your 64-bit perform method
	// 4: flags to alter how the signal chain handles your object -- just pass 0
	// 5: a generic pointer that you can use to pass any additional data to your perform method
	
	object_method(dsp64, gensym("dsp_add64"), x, specRolloff_perform64, 0, NULL);
    
    
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
        
		// freqs for each bin based on current window size and sample rate
		for(i=0; i<x->window; i++)
			x->binFreqs[i] = (x->sr/x->window)*i;
        
    	post("specRolloff~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)x->sr, (int)x->n);
	};
}


//***********************************************************************************************

// this is the 32-bit perform method for Max 5 and earlier
t_int *specRolloff_perform(t_int *w)
{
    int i, n;
    double timef;
	// DO NOT CALL post IN HERE, but you can call defer_low (not defer)
	
	// args are in a vector, sized as specified in specRolloff_dsp method
	// w[0] contains &specRolloff_perform, so we start at w[1]
	t_specRolloff *x = (t_specRolloff *)(w[1]);
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
void specRolloff_perform64(t_specRolloff *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
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

