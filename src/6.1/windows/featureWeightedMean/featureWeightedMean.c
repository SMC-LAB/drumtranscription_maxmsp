#include "ext.h"							// standard Max include, always required
#include "ext_obex.h"						// required for new style Max object
#include "m_memory.h"	
#define MAXCLUSTERMEMS 8192
#define t_float float

////////////////////////// object struct
typedef struct instance
{
    float *instance;
} t_instance;

typedef struct _featureWeightedMean
{
	t_object  x_obj;
    t_atom *x_listOut;
    t_instance *instances;
    int spew;
    int featureLength;
    int currentFrame;
    int numInstances;
    int numFrames;
    int weight;
    
    void *featureList;
} t_featureWeightedMean;


///////////////////////// function prototypes
static void featureWeightedMean_accum(t_featureWeightedMean *x, t_symbol *s, int argc, t_atom *argv);
static void featureWeightedMean_numFrames(t_featureWeightedMean *x, int num);
static void featureWeightedMean_length(t_featureWeightedMean *x,int len);
static void featureWeightedMean_print(t_featureWeightedMean *x);
static void featureWeightedMean_clear(t_featureWeightedMean *x);
static void featureWeightedMean_bang(t_featureWeightedMean  *x);
static void featureWeightedMean_free(t_featureWeightedMean *x);
static void *featureWeightedMean_new(t_symbol *s, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *featureWeightedMean_class;


int C74_EXPORT main(void)
{	
	t_class *c;
	
	c = class_new("featureWeightedMean", (method)featureWeightedMean_new, (method)featureWeightedMean_free, (long)sizeof(t_featureWeightedMean),0L /* leave NULL!! */,  A_GIMME, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	featureWeightedMean_class = c;
    
    addbang((method)featureWeightedMean_bang);
    class_addmethod(featureWeightedMean_class, (method)featureWeightedMean_bang, "bang", 0);
    
    class_addmethod(featureWeightedMean_class, (method)featureWeightedMean_accum, "list", A_GIMME, 0);
    
    
    
    class_addmethod(featureWeightedMean_class, 
                    (method)featureWeightedMean_clear,
                    "clear",
                    0
                    );
    
    class_addmethod(featureWeightedMean_class, 
                    (method)featureWeightedMean_numFrames,
                    "num_frames",
                    A_LONG,
                    0
                    );
    
    
    class_addmethod(featureWeightedMean_class, 
                    (method)featureWeightedMean_length,
                    "length",
                    A_LONG,
                    0
                    );
        
    class_addmethod(featureWeightedMean_class, 
                    (method)featureWeightedMean_bang,
                    "bang",
                    0
                    );  
    
    class_addmethod(featureWeightedMean_class, 
                    (method)featureWeightedMean_print,
                    "print",
                    0
                    );
    

	post("featureWeightedMean version 0.1.0");
	return 0;
}

static void featureWeightedMean_accum(t_featureWeightedMean *x, t_symbol *s, int argc, t_atom *argv)
{
	int i, j;
    float sum,w;
    
	if(x->featureLength != argc)
		error("featureWeightedMean: input length does not match current length setting. input ignored.");
	else
		for(i=0; i<x->featureLength; i++)
			x->instances[x->currentFrame].instance[i] = atom_getfloat(argv+i);
    
	x->currentFrame++;
    
	if (x->currentFrame==x->numFrames) 
	{
        for(i=0; i<(x->featureLength-1); i++)
        {
            sum = 0; 
            w = 0;
            for(j=0; j<x->currentFrame; j++) 
            {
                sum = sum + x->instances[j].instance[i] * x->instances[j].instance[x->featureLength-1];
                w = w + x->instances[j].instance[x->featureLength-1];
            }
            sum = sum / w;
            atom_setfloat(x->x_listOut+i, sum); 
        }
               
		outlet_list(x->featureList, 0, x->featureLength-1, x->x_listOut);
        
		x->currentFrame = (x->currentFrame==x->numFrames)?0:x->currentFrame;
	}
}

static void featureWeightedMean_bang(t_featureWeightedMean *x)
{
    int i, j;
    float sum,w;
    
    if (x->currentFrame>0) 
    {
        for(i=0; i<(x->featureLength-1); i++)
        {
            sum = 0; 
            w = 0;
            for(j=0; j<x->currentFrame; j++) 
            {
                sum = sum + x->instances[j].instance[i] * x->instances[j].instance[x->featureLength-1];
                w = w + x->instances[j].instance[x->featureLength-1];
            }
            sum = sum / w;
            atom_setfloat(x->x_listOut+i, sum); 
        }
        
        outlet_list(x->featureList, 0, x->featureLength-1, x->x_listOut);
        
        x->currentFrame = 0;
    
    }
}

static void featureWeightedMean_clear(t_featureWeightedMean *x)
{
	int i, j;
    
	// free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
	t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
    
	x->currentFrame = 0;
    
    for(i=0; i<(x->featureLength-1); i++)
        atom_setfloat(x->x_listOut+i, 0.0);
    
    x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
    
	for(i=0; i<x->numFrames; i++)
		x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
    
	for(i=0; i<x->numFrames; i++)
		for(j=0; j<x->featureLength; j++)
			x->instances[i].instance[j] = 0.0;
}

static void featureWeightedMean_numFrames(t_featureWeightedMean *x, int num)
{
	int i, j;
    
	if(num)
	{
        //x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, (x->featureLength-1)*sizeof(t_atom), (x->featureLength*num)*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->currentFrame = 0;
		x->numFrames = num;
        
        for(i=0; i<(x->featureLength-1); i++)
	        atom_setfloat(x->x_listOut+i, 0.0);
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureWeightedMean_length(t_featureWeightedMean *x, int len)
{
	int i, j;
    
	if(len)
	{
        x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, (x->featureLength-1)*sizeof(t_atom), (len*x->numFrames)*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		x->featureLength = len;
		x->currentFrame = 0;
        
        for(i=0; i<(x->featureLength-1); i++)
	        atom_setfloat(x->x_listOut+i, 0.0);
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureWeightedMean_print(t_featureWeightedMean *x)
{
    post("averaging with weights %i vectors with %i features", x->numFrames, x->featureLength);
}

static void featureWeightedMean_free(t_featureWeightedMean *x)
{
	int i;      
    
    // free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
    t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
    
    // free listOut memory
	t_freebytes_(x->x_listOut, (x->featureLength-1)*sizeof(t_atom));
    
}

static void *featureWeightedMean_new(t_symbol *s, long argc, t_atom *argv)
{
    t_featureWeightedMean *x = NULL;
    int i,j;
    t_atom *ap;

	if ((x = (t_featureWeightedMean *)object_alloc(featureWeightedMean_class))) {
        
        x->featureLength = 2;
        x->numFrames = 2;
        x->currentFrame = 0;
        // increment ap each time to get to the next atom
        for (i = 0, ap = argv; i < argc; i++, ap++) 
        {
            if (atom_gettype(ap) == A_LONG) 
            {
                if (i==0) x->featureLength = atom_getlong(ap);
                if (i==1) x->numFrames = atom_getlong(ap);
            }
            else post("%ld: unknown atom type (%ld)", i+1, atom_gettype(ap));
        }
        
        
        x->featureList = listout(x);
     
        x->x_listOut = (t_atom *)t_getbytes_((x->featureLength-1)*sizeof(t_atom));
        x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
        for(i=0; i<(x->featureLength-1); i++)
            atom_setfloat(x->x_listOut+i, 0.0);
        
        for(i=0; i<x->numFrames; i++)
            x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
        for(i=0; i<x->numFrames; i++)
            for(j=0; j<x->featureLength; j++)
                x->instances[i].instance[j] = 0.0;

    }
    
    return (x);
}

