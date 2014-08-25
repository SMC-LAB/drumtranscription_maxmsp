/**
 THIS CODE WAS ADAPTED FOR MAXMSP FROM THE AUBIO LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
 @file
 aubioOnset~: a simple audio object for Max
 @ingroup examples	
 */

#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects
#include "aubio.h"			// required for Aubio
#include "m_memory.h"
#define t_float float

// struct to represent the object's state
typedef struct _aubioOnset {
	t_pxobject		ob;			// the object itself (t_pxobject in MSP instead of t_object)
    t_float sr;
    t_float n;
    
    t_float threshold;
    t_float silence;
    int minioi;  
    int pos;                    /*frames%dspblocksize */
    int bufsize;
    int hopsize;
    char* detectionFunction;
    
    fvec_t *in;
    fvec_t *out;
    aubio_onset_t *o;
    
    void *onsetbang;
    
} t_aubioOnset;


// method prototypes
static void aubioOnset_threshold (t_aubioOnset * x, double f);
static void aubioOnset_silence (t_aubioOnset * x, int s);
static void aubioOnset_minioi (t_aubioOnset * x, int m);
static void aubioOnset_debug (t_aubioOnset * x);

void *aubioOnset_new(t_symbol *s, long argc, t_atom *argv);
void aubioOnset_free(t_aubioOnset *x);
void aubioOnset_dsp(t_aubioOnset *x, t_signal **sp, short *count);
void aubioOnset_dsp64(t_aubioOnset *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
t_int *aubioOnset_perform(t_int *w);
void aubioOnset_perform64(t_aubioOnset *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *aubioOnset_class = NULL;


//***********************************************************************************************

int C74_EXPORT main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.
    
	t_class *c = class_new("aubioOnset~", (method)aubioOnset_new, (method)dsp_free, (long)sizeof(t_aubioOnset), 0L, A_GIMME, 0);
    
    class_dspinit(c);
	class_register(CLASS_BOX, c);
	aubioOnset_class = c;
	
	class_addmethod(c, (method)aubioOnset_dsp,		"dsp",		A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 5 and earlier
	class_addmethod(c, (method)aubioOnset_dsp64,		"dsp64",	A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
	//class_addmethod(c, (method)aubioOnset_assist,	"assist",	A_CANT, 0);
    
    class_addmethod(aubioOnset_class, 
                    (method)aubioOnset_threshold,
                    "threshold",
                    A_DEFFLOAT,
                    0
                    );	
    
    class_addmethod(aubioOnset_class, 
                    (method)aubioOnset_silence,
                    "silence",
                    A_DEFLONG,
                    0
                    );	
    
    class_addmethod(aubioOnset_class, 
                    (method)aubioOnset_minioi,
                    "minioi",
                    A_DEFLONG,
                    0
                    );	 
    
    class_addmethod(aubioOnset_class, 
                    (method)aubioOnset_debug,
                    "debug",
                    0
                    );
    
	return 0;
}


void *aubioOnset_new(t_symbol *s, long argc, t_atom *argv)
{
	t_aubioOnset *x = (t_aubioOnset *)object_alloc(aubioOnset_class);
    t_atom *ap;
    int i, isPow2, argcount=0;
	//s=s;
    
	if (x) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
        // use 0 if you don't need inlets
		//outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)        
        
        x->onsetbang = bangout(x);
        
        x->sr = 44100.0;
        x->n = 64.0;
        
        x->threshold = 0.3;
        x->silence = -70;
        x->minioi = 4;
        x->bufsize = 1024;
        x->hopsize = 512;
        x->detectionFunction = "complex";
        
        // increment ap each time to get to the next atom
        for (i = 0, ap = argv; i < argc; i++, ap++) 
        {
            switch (atom_gettype(ap)) {
                case A_LONG:
                    if (atom_getlong(ap)<0)
                    {
                        post("%ld: silence threshold %ld",i+1,atom_getlong(ap));
                        x->silence = (atom_getlong(ap) < -120) ? -120 : (atom_getlong(ap) > 0) ? 0 : atom_getlong(ap);
                    }
                    else 
                    {
                        if (argcount == 0) {
                            post("%ld: bufsize %ld",i+1,atom_getlong(ap));
                            x->bufsize = atom_getlong(ap);
                            argcount = argcount + 1;       
                        }
                        else {
                            post("%ld: hopsize %ld",i+1,atom_getlong(ap));
                            x->hopsize = atom_getlong(ap);
                        }
                    }
                    break;
                case A_FLOAT:
                    post("%ld: threshold %.2f",i+1,atom_getfloat(ap));
                    x->threshold = (atom_getfloat(ap) < 1e-5) ? 0.1 : (atom_getfloat(ap) > 0.999) ? 0.999 : atom_getfloat(ap);
                    break;
                case A_SYM:
                    post("%ld: onset detection function %s",i+1, atom_getsym(ap)->s_name);
                    x->detectionFunction = atom_getsym(ap)->s_name;
                    break;
                default:
                    post("%ld: unknown atom type (%ld)", i+1, atom_gettype(ap));
                    break;
            }
        }
        
        isPow2 = (int)x->bufsize && !( ((int)x->bufsize-1) & (int)x->bufsize );            
        if(!isPow2)
        {
            error("requested buffer size is not a power of 2. default value of 1024 used instead");
            x->bufsize = 1024;
        }
        isPow2 = (int)x->hopsize && !( ((int)x->hopsize-1) & (int)x->hopsize );            
        if(!isPow2)
        {
            error("requested hop size is not a power of 2. default value of 1024 used instead");
            x->hopsize = x->bufsize / 4;
        }
        
        
        if (strcmp(x->detectionFunction,"hfc") == 0) x->o=new_aubio_onset(aubio_onset_hfc,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"energy") == 0) x->o=new_aubio_onset(aubio_onset_energy,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"phase") == 0) x->o=new_aubio_onset(aubio_onset_phase,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"complex") == 0) x->o=new_aubio_onset(aubio_onset_complex,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"specdiff") == 0) x->o=new_aubio_onset(aubio_onset_specdiff,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"kl") == 0) x->o=new_aubio_onset(aubio_onset_kl,x->bufsize, x->hopsize, 1);
        else if (strcmp(x->detectionFunction,"mkl") == 0) x->o=new_aubio_onset(aubio_onset_mkl,x->bufsize, x->hopsize, 1);
        else x->o=new_aubio_onset(aubio_onset_complex,x->bufsize, x->hopsize, 1);             
        
        
        x->in = (fvec_t *) new_fvec (x->hopsize, 1);
        x->out = (fvec_t *) new_fvec (1, 1);
        
        aubio_onset_set_threshold(x->o,x->threshold); 
        aubio_onset_set_silence(x->o,x->silence);  
        aubio_onset_set_minioi(x->o,x->minioi);  
            
        post("aubioOnset~: version 0.3");
        
	}
	return (x);
}


// NOT CALLED!, we use dsp_free for a generic free function
void aubioOnset_free(t_aubioOnset *x) 
{
	// free method
    

}
//***********************************************************************************************


static void aubioOnset_threshold (t_aubioOnset * x, double f)
{
    //x->threshold = (f < 1e-5) ? 0.1 : (f > 0.999) ? 0.999 : f;
    //aubio_onset_set_threshold(x->o,f);
    if (f < 1e-5)
        x->threshold = 0.1;
    else if (f > 0.999)
        x->threshold = 0.9;
    else
        x->threshold = f;
    aubio_onset_set_threshold(x->o, x->threshold);
    post ("aubioonsethfc~ threshold:\t%f", x->threshold);
}

static void aubioOnset_silence (t_aubioOnset * x, int s)
{
    x->silence = (s < -120) ? -120 : (s > 0.) ? 0. : s;
    aubio_onset_set_silence(x->o,x->silence); 
    post ("aubioonsethfc~ silence threshold:\t%f", x->silence);
}

static void aubioOnset_minioi (t_aubioOnset * x, int m)
{
    x->minioi = (m < 4) ? 4 : (int)m;
    aubio_onset_set_minioi(x->o,x->minioi); 
    post ("aubioonsethfc~ minioi:\t%d", x->minioi);
}

static void aubioOnset_debug (t_aubioOnset * x)
{
    post ("aubioonsethfc~ bufsize:\t%d", x->bufsize);
    post ("aubioonsethfc~ hopsize:\t%d", x->hopsize);
    post ("aubioonsethfc~ threshold:\t%f", x->threshold);
    post ("aubioonsethfc~ silence threshold:\t%f", x->silence);
    post ("aubioonsethfc~ minioi:\t%d", x->minioi);
    post ("aubioonsethfc~ audio in:\t%f", x->in->data[0][0]);
    post ("aubioonsethfc~ onset:\t%f", x->out->data[0][0]);
}

//***********************************************************************************************
// this function is called when the DAC is enabled, and "registers" a function for the signal chain in Max 5 and earlier.
// In this case we register the 32-bit, "aubioOnset_perform" method.
void aubioOnset_dsp(t_aubioOnset *x, t_signal **sp, short *count)
{
	//post("my sample rate is: %f", sp[0]->s_sr);
	
	// dsp_add
	// 1: (t_perfroutine p) perform method
	// 2: (long argc) number of args to your perform method
	// 3...: argc additional arguments, all must be sizeof(pointer) or long
	// these can be whatever, so you might want to include your object pointer in there
	// so that you have access to the info, if you need it.
	dsp_add(aubioOnset_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void aubioOnset_dsp64(t_aubioOnset *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	//post("my sample rate is: %f", samplerate);
	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the dsp_add64 arguments are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: a pointer to your object
	// 3: a pointer to your 64-bit perform method
	// 4: flags to alter how the signal chain handles your object -- just pass 0
	// 5: a generic pointer that you can use to pass any additional data to your perform method
	
	object_method(dsp64, gensym("dsp_add64"), x, aubioOnset_perform64, 0, NULL);
    
}


//***********************************************************************************************

// this is the 32-bit perform method for Max 5 and earlier
t_int *aubioOnset_perform(t_int *w)
{
    int j, n;
	// DO NOT CALL post IN HERE, but you can call defer_low (not defer)
	
	// args are in a vector, sized as specified in aubioOnset_dsp method
	// w[0] contains &aubioOnset_perform, so we start at w[1]
	t_aubioOnset *x = (t_aubioOnset *)(w[1]);
	//t_float *inL = (t_float *)(w[2]);
	//t_float *outL = (t_float *)(w[3]);
    t_sample *in = (t_float *)(w[2]);
	n = (int)w[3];
    
    for (j = 0; j < n; j++) {
        /* write input to datanew */
        fvec_write_sample (x->in, in[j], 0, x->pos);
        /*time to do something */
        if (x->pos == x->hopsize - 1) {
            /* block loop */
            aubio_onset (x->o, x->in, x->out);
            if (fvec_read_sample (x->out, 0, 0) > 0.) {
                outlet_bang (x->onsetbang);
            }
            /* end of block loop */
            x->pos = -1;              /* so it will be zero next j loop */
        }
        x->pos++;
    }
    
	// you have to return the NEXT pointer in the array OR MAX WILL CRASH
	return w + 4;
}


// this is 64-bit perform method for Max 6
void aubioOnset_perform64(t_aubioOnset *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    int j, n;
	//t_double *inL = ins[0];		// we get audio for each inlet of the object from the **ins argument
    t_sample *in = ins[0];
	//t_double *outL = outs[0];	// we get audio for each outlet of the object from the **outs argument
	n = sampleframes;

    for (j = 0; j < n; j++) {
        /* write input to datanew */
        fvec_write_sample (x->in, in[j], 0, x->pos);
        /*time to do something */
        if (x->pos == x->hopsize - 1) {
            /* block loop */
            aubio_onset (x->o, x->in, x->out);
            if (fvec_read_sample (x->out, 0, 0) > 0.) {
                outlet_bang (x->onsetbang);
            }
            /* end of block loop */
            x->pos = -1;              /* so it will be zero next j loop */
        }
        x->pos++;
    }

	
	// this perform method simply copies the input to the output, offsetting the value
	//while (n--)
	//	*outL++ = *inL++ + x->offset;
}

