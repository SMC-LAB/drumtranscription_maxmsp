/*
  THIS CODE WAS ADAPTED FOR MAXMSP FROM THE TIMBREID PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
 
 barkSpec~ - A Bark Frequency Spectrum Analysis external.
 
 Copyright 2009 William Brent
 
 This file is part of timbreID.
 
 timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
 timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 
 version 0.2.9D, November 28, 2011
 
 • 0.2.9D finally got memory resizing in tIDLib header to work. Also including Nyquist now just for good measure.  Bringing back overlap, since it's the only way to inform the object of something like [block~ 1024 4], where the sampling rate will appear as 4x its actual rate.  And sampling rate affects filterbank construction....
 • 0.2.9C lost track of changes, but bottom line is that x->x_filterbankSizes memory resizing isn't working in create_filterbank.
 • 0.2.9B changing variable names to cut down on underscores
 • 0.2.9 employing the timbreID.h header in preparation for timbreID-0.6.0.  barkSpec~.c is the first attempt to integrate this.  Going to get rid of the x->overlap variable too.  Not necessary anymore.
 • 0.2.8 fixed a memory leak: wasn't freeing memory for x_filterbankSizes. made power spectrum computation the default, and changed the squaring function to a magnitude function instead.  in the case that power spectrum is used, this saves needless computation of sqrt and subsequent squaring. wherever possible, using t_getbytes() directly instead of getting 0 bytes and resizing. corrected createFilterbank declaration in _setup to show that it only has one argument, not two.
 • 0.2.7 fixed a memory issue in _make_filterbank(): was freeing window*sizeof(float) bytes instead of window_half*sizeof(float) bytes
 • 0.2.6 decided to go with a switch for either summing power in filters, or averaging. default is summing.
 • 0.2.5 changed filterbank_multiply so that sum of power in each filter is divided by the number of points in the filter.
 • 0.2.4 added an ifndef M_PI for guaranteed windows compilation
 • 0.2.3 adds a #define M_PI for windows compilation, and declares all functions except _setup static
 • 0.2.2 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).
 • 0.2.1 uses dynamic memory allocation for filter widths
 • 0.2.0 implements mayer_realfft
 • 0.1.9 added normalization option
 
 */

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include "tIDLib.h"
#define t_float float

// struct to represent the object's state
typedef struct _barkSpec {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
	//double			offset; 	// the value of a property of our object
    t_float sr;
    t_float n;
	int windowFunction;
    int overlap;
    int powerSpectrum;
    int window;
    int normalize;
	int filterAvg;
    int sizeFilterFreqs;
	int numFilters;
    double lastDspTime;
    t_sample *signal_R;
    t_float barkSpacing;
	t_float *x_filterFreqs;
	t_filter *x_filterbank;
    t_float *blackman;
    t_float *cosine;
    t_float *hamming;
    t_float *hann;    
    t_float x_f;
    
    void *x_featureList;

} t_barkSpec;


// method prototypes
static void barkSpec_tilde_bang(t_barkSpec *x);
static void barkSpec_tilde_window(t_barkSpec *x, int w);
static void barkSpec_tilde_overlap(t_barkSpec *x, int o);
static void barkSpec_tilde_windowFunction(t_barkSpec *x, int f);
static void barkSpec_tilde_powerSpectrum(t_barkSpec *x, int spec);
static void barkSpec_tilde_createFilterbank(t_barkSpec *x, t_float bs);
static void barkSpec_tilde_filterFreqs(t_barkSpec *x);
static void barkSpec_tilde_print(t_barkSpec *x);
static void barkSpec_tilde_hat(t_barkSpec *x, int filt);
static void barkSpec_tilde_normalize(t_barkSpec *x, int norm);
static void barkSpec_tilde_filterAvg(t_barkSpec *x, int avg);

void *barkSpec_new(t_symbol *s, long argc, t_atom *argv);
void barkSpec_free(t_barkSpec *x);
void barkSpec_dsp(t_barkSpec *x, t_signal **sp, short *count);
void barkSpec_dsp64(t_barkSpec *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *barkSpec_perform(t_int *w);
void barkSpec_perform64(t_barkSpec *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *barkSpec_class = NULL;


//***********************************************************************************************

int C74_EXPORT main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	t_class *c = class_new("barkSpec~", (method)barkSpec_new, (method)dsp_free, (long)sizeof(t_barkSpec), 0L, A_GIMME, 0);
    
    class_dspinit(c);
	class_register(CLASS_BOX, c);
	barkSpec_class = c;
	
	class_addmethod(c, (method)barkSpec_dsp,		"dsp",		A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 5 and earlier
	class_addmethod(c, (method)barkSpec_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
	//class_addmethod(c, (method)barkSpec_assist,	"assist",	A_CANT, 0);
    
    
	addbang((method)barkSpec_tilde_bang);
    class_addmethod(barkSpec_class, 
                    (method)barkSpec_tilde_bang,
                    "bang",
                    0
                    ); 
    
    class_addmethod(barkSpec_class, 
                    (method)barkSpec_tilde_filterFreqs,
                    "filter_freqs",
                    0
                    ); 
    
    class_addmethod(barkSpec_class, 
                    (method)barkSpec_tilde_print,
                    "print",
                    0
                    );     
    
	class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_window,
                    "window",
                    A_LONG,
                    0
                    );
    
	class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_overlap,
                    "overlap",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_windowFunction,
                    "window_function",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_powerSpectrum,
                    "power_spectrum",
                    A_DEFLONG,
                    0
                    );
    
    class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_createFilterbank,
                    "filterbank",
                    A_DEFFLOAT,
                    0
                    );
    
	class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_hat,
                    "hat",
                    A_DEFFLOAT,
                    0
                    );
    
    class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_filterAvg,
                    "filter_avg",
                     A_DEFLONG,
                    0
                    );
    
    class_addmethod(
                    barkSpec_class,
                    (method)barkSpec_tilde_normalize,
                    "normalize",
                     A_DEFLONG,
                    0
                    );
	


	return 0;
}


void *barkSpec_new(t_symbol *s, long argc, t_atom *argv)
{
	t_barkSpec *x = (t_barkSpec *)object_alloc(barkSpec_class);
    int i, isPow2;
    double timef = 0.0;
	//s=s;
    
	if (x) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
										// use 0 if you don't need inlets
		//outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
		//x->offset = 0.0;
        
        
        x->x_featureList = listout(x);
        //x->x_featureList = outlet_new(&x->x_obj, &s_float);        
        
        if(argc > 1)
        {
            x->window = atom_getlong(argv);  // should perform a check for >64 && power of two
            isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
           
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead.");
                x->window = 1024;
            };
            
            x->barkSpacing = atom_getfloat(argv+1);
        }
        else if(argc > 0)
        {
            x->window = atom_getlong(argv);
            isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
            
            if(!isPow2)
            {
                error("requested window size is not a power of 2. default value of 1024 used instead.");
                x->window = 1024;
            };
            
            x->barkSpacing = 0.5;
        }
        else
        {
            x->window = 1024;
            x->barkSpacing = 0.5;
        }
        
        x->sr = 44100.0;
        x->n = 64.0;
        x->overlap = 1;
        x->windowFunction = 4; // 4 is hann window
        x->powerSpectrum = 0; // choose mag (0) or power (1) spec
        x->normalize = 1; // this is generally a good thing, but should be off for concatenative synth
        x->filterAvg = 0;
        //x->lastDspTime = clock_getlogicaltime();
        clock_getftime(&timef);
        x->lastDspTime = timef;
        x->sizeFilterFreqs = 0;
        x->numFilters = 0; // this is just an init size that will be updated in createFilterbank anyway.
        
        x->signal_R = (t_sample *)t_getbytes_((x->window+x->n)*sizeof(t_sample));
        
        // initialize signal buffer
        for(i=0; i<x->window+x->n; i++)
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
        
        // grab memory
        x->x_filterbank = (t_filter *)t_getbytes_(0);
        x->x_filterFreqs = (t_float *)t_getbytes_(0);
        
        x->numFilters = tIDLib_getBarkFilterSpacing(&x->x_filterFreqs, x->sizeFilterFreqs, x->barkSpacing, x->sr);
        
        x->sizeFilterFreqs = x->numFilters+2;
        
        tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, 0, x->numFilters, x->window, x->sr);
        
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void barkSpec_free(t_barkSpec *x) 
{
    int i;
    
	// free the input buffer memory
    t_freebytes_(x->signal_R, (x->window+x->n)*sizeof(t_sample));
    
    
    // free filterFreqs memory
	t_freebytes_(x->x_filterFreqs, x->sizeFilterFreqs*sizeof(t_float));
    
	// free the filterbank memory
	for(i=0; i<x->numFilters; i++)
		t_freebytes_(x->x_filterbank[i].filter, x->x_filterbank[i].size*sizeof(t_float));
    
    t_freebytes_(x->x_filterbank, x->numFilters*sizeof(t_filter));
    
	// free the window memory
	t_freebytes_(x->blackman, x->window*sizeof(t_float));
	t_freebytes_(x->cosine, x->window*sizeof(t_float));
	t_freebytes_(x->hamming, x->window*sizeof(t_float));
	t_freebytes_(x->hann, x->window*sizeof(t_float));
}


//***********************************************************************************************

static void barkSpec_tilde_bang(t_barkSpec *x)
{
    int i, j, window, windowHalf, bangSample;
	t_atom *listOut;
	t_sample *signal_R, *signal_I;
	t_float *windowFuncPtr;
	double currentTime, timef = 0.0;
    
	window = x->window;
	windowHalf = window*0.5;
    
	// create local memory
	listOut = (t_atom *)t_getbytes_(x->numFilters*sizeof(t_atom));
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
    
	// power spectrum sometimes generates lower scores than magnitude. make it optional.
	if(!x->powerSpectrum)
		tIDLib_mag(windowHalf+1, signal_R);
    
	tIDLib_filterbankMultiply(signal_R, x->normalize, x->filterAvg, x->x_filterbank, x->numFilters);
    
	for(i=0; i<x->numFilters; i++)
		atom_setfloat(listOut+i, signal_R[i]);
    
	outlet_list(x->x_featureList, 0, x->numFilters, listOut);
    
	// free local memory
	t_freebytes_(listOut, x->numFilters*sizeof(t_atom));
	t_freebytes_(signal_R, window*sizeof(t_sample));
	t_freebytes_(signal_I, (windowHalf+1)*sizeof(t_sample));    
    
}

static void barkSpec_tilde_createFilterbank(t_barkSpec *x, t_float bs)
{
	int oldNumFilters;
    
	x->barkSpacing = bs;
	oldNumFilters = x->numFilters;
    
	x->numFilters = tIDLib_getBarkFilterSpacing(&x->x_filterFreqs, x->sizeFilterFreqs, x->barkSpacing, x->sr);
    
	x->sizeFilterFreqs = x->numFilters+2;
    
	tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->numFilters, x->window, x->sr);
}


static void barkSpec_tilde_filterFreqs(t_barkSpec *x)
{
	int i;
    
	for(i=0; i<x->numFilters+2; i++)
		post("filterFreqs[%i]: %f", i, x->x_filterFreqs[i]);
}


static void barkSpec_tilde_print(t_barkSpec *x)
{
	post("samplerate: %f", x->sr);
	post("window: %i", x->window);
    
	post("bark spacing: %f", x->barkSpacing);
	post("power spectrum: %i", x->powerSpectrum);
	post("normalize: %i", x->normalize);
	post("filter averaging: %i", x->filterAvg);
	post("window function: %i", x->windowFunction);
    
	post("no. of filters: %i", x->numFilters);
}

static void barkSpec_tilde_hat(t_barkSpec *x, int filt)
{
	int i, idx;
    
	idx = filt;
    
	if(idx>=x->numFilters)
		error("filter %i does not exist.", idx);
	else if(idx < 0)
		error("filter %i does not exist.", idx);
	else
	{
		post("size[%i]: %i", idx, x->x_filterbank[idx].size);
		for(i=0; i<x->x_filterbank[idx].size; i++)
			post("val %i: %f", i, x->x_filterbank[idx].filter[i]);
        
		post("idxLo: %i, idxHi: %i", x->x_filterbank[idx].indices[0], x->x_filterbank[idx].indices[1]);
	}
}



static void barkSpec_tilde_window(t_barkSpec *x, int w)
{
    int i, window, isPow2, oldNumFilters;
	
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
        
		x->window = (t_float)window;
        
		// re-init window functions
		tIDLib_blackmanWindow(x->blackman, x->window);
		tIDLib_cosineWindow(x->cosine, x->window);
		tIDLib_hammingWindow(x->hamming, x->window);
		tIDLib_hannWindow(x->hann, x->window);
        
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
        
        oldNumFilters = x->numFilters;
		x->numFilters = tIDLib_getBarkFilterSpacing(&x->x_filterFreqs, x->sizeFilterFreqs, x->barkSpacing, x->sr);
        
		x->sizeFilterFreqs = x->numFilters+2;
        
		tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->numFilters, x->window, x->sr);

		
		post("window size: %i", (int)x->window);
	}

}


static void barkSpec_tilde_overlap(t_barkSpec *x, int o)
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


static void barkSpec_tilde_windowFunction(t_barkSpec *x, int f)
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

static void barkSpec_tilde_filterAvg(t_barkSpec *x, int avg)
{
    avg = (avg<0)?0:avg;
    avg = (avg>1)?1:avg;
	x->filterAvg = avg;
    
	if(x->filterAvg)
		post("averaging energy in filters.");
	else
		post("summing energy in filters.");
}


// magnitude spectrum == 0, power spectrum == 1
static void barkSpec_tilde_powerSpectrum(t_barkSpec *x, int spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->powerSpectrum = spec;
    
	if(x->powerSpectrum)
		post("using power spectrum");
	else
		post("using magnitude spectrum");
}

static void barkSpec_tilde_normalize(t_barkSpec *x, int norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
	x->normalize = norm;
    
	if(x->normalize)
		post("spectrum normalization ON.");
	else
		post("spectrum normalization OFF.");
}



//***********************************************************************************************

// this function is called when the DAC is enabled, and "registers" a function for the signal chain in Max 5 and earlier.
// In this case we register the 32-bit, "barkSpec_perform" method.
void barkSpec_dsp(t_barkSpec *x, t_signal **sp, short *count)
{
    int i, oldNumFilters;
    double timef;
	//post("my sample rate is: %f", sp[0]->s_sr);
	
	// dsp_add
	// 1: (t_perfroutine p) perform method
	// 2: (long argc) number of args to your perform method
	// 3...: argc additional arguments, all must be sizeof(pointer) or long
	// these can be whatever, so you might want to include your object pointer in there
	// so that you have access to the info, if you need it.
	dsp_add(barkSpec_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
    
    
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
                
    	oldNumFilters = x->numFilters;
		x->numFilters = tIDLib_getBarkFilterSpacing(&x->x_filterFreqs, x->sizeFilterFreqs, x->barkSpacing, x->sr);
        
		x->sizeFilterFreqs = x->numFilters+2;
        
		tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->numFilters, x->window, x->sr);
        
		post("barkSpec~: window size: %i. sampling rate: %i, block size: %i", (int)x->window, (int)x->sr, (int)x->n);
	};
}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void barkSpec_dsp64(t_barkSpec *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
    int i, oldNumFilters;
    double timef;
	//post("my sample rate is: %f", samplerate);
	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the dsp_add64 arguments are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: a pointer to your object
	// 3: a pointer to your 64-bit perform method
	// 4: flags to alter how the signal chain handles your object -- just pass 0
	// 5: a generic pointer that you can use to pass any additional data to your perform method
	
	object_method(dsp64, gensym("dsp_add64"), x, barkSpec_perform64, 0, NULL);
    
    
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
        
    	oldNumFilters = x->numFilters;
		x->numFilters = tIDLib_getBarkFilterSpacing(&x->x_filterFreqs, x->sizeFilterFreqs, x->barkSpacing, x->sr);
        
		x->sizeFilterFreqs = x->numFilters+2;
        
		tIDLib_createFilterbank(x->x_filterFreqs, &x->x_filterbank, oldNumFilters, x->numFilters, x->window, x->sr);
        
		post("barkSpec~: window size: %i. sampling rate: %i, block size: %i", (int)x->window, (int)x->sr, (int)x->n);
	};
}


//***********************************************************************************************

// this is the 32-bit perform method for Max 5 and earlier
t_int *barkSpec_perform(t_int *w)
{
    int i, n;
    double timef;
	// DO NOT CALL post IN HERE, but you can call defer_low (not defer)
	
	// args are in a vector, sized as specified in barkSpec_dsp method
	// w[0] contains &barkSpec_perform, so we start at w[1]
	t_barkSpec *x = (t_barkSpec *)(w[1]);
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
void barkSpec_perform64(t_barkSpec *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
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

