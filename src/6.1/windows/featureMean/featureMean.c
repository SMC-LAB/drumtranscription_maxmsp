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

typedef struct _featureMean
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
} t_featureMean;


///////////////////////// function prototypes
static void featureMean_accum(t_featureMean *x, t_symbol *s, int argc, t_atom *argv);
static void featureMean_numFrames(t_featureMean *x, int num);
static void featureMean_length(t_featureMean *x,int len);
static void featureMean_print(t_featureMean *x);
static void featureMean_clear(t_featureMean *x);
static void featureMean_bang(t_featureMean  *x);
static void featureMean_free(t_featureMean *x);
static void *featureMean_new(t_symbol *s, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *featureMean_class;


int C74_EXPORT main(void)
{	
	t_class *c;
	
	c = class_new("featureMean", (method)featureMean_new, (method)featureMean_free, (long)sizeof(t_featureMean),0L /* leave NULL!! */,  A_GIMME, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	featureMean_class = c;
    
    addbang((method)featureMean_bang);
    class_addmethod(featureMean_class, (method)featureMean_bang, "bang", 0);
    
    class_addmethod(featureMean_class, (method)featureMean_accum, "list", A_GIMME, 0);
    
    
    
    class_addmethod(featureMean_class, 
                    (method)featureMean_clear,
                    "clear",
                    0
                    );
    
    class_addmethod(featureMean_class, 
                    (method)featureMean_numFrames,
                    "num_frames",
                    A_LONG,
                    0
                    );
    
    
    class_addmethod(featureMean_class, 
                    (method)featureMean_length,
                    "length",
                    A_LONG,
                    0
                    );
        
    class_addmethod(featureMean_class, 
                    (method)featureMean_bang,
                    "bang",
                    0
                    );  
    
    class_addmethod(featureMean_class, 
                    (method)featureMean_print,
                    "print",
                    0
                    );
    

	post("featureMean version 0.1.0");
	return 0;
}

static void featureMean_accum(t_featureMean *x, t_symbol *s, int argc, t_atom *argv)
{
	int i, j;
    float sum;
    
	if(x->featureLength != argc)
		error("featureMean: input length does not match current length setting. input ignored.");
	else
		for(i=0; i<x->featureLength; i++)
			x->instances[x->currentFrame].instance[i] = atom_getfloat(argv+i);
    
	x->currentFrame++;
    
	if (x->currentFrame==x->numFrames) 
	{
        for(i=0; i<x->featureLength; i++)
        {
            sum = 0;
            for(j=0; j<x->currentFrame; j++) 
            {
                sum = sum + x->instances[j].instance[i];
            }
            sum = sum / x->currentFrame;
            atom_setfloat(x->x_listOut+i, sum); 
        }
               
		outlet_list(x->featureList, 0, x->featureLength, x->x_listOut);
        
		x->currentFrame = (x->currentFrame==x->numFrames)?0:x->currentFrame;
	}
}

static void featureMean_bang(t_featureMean *x)
{
    int i, j;
    float sum;
    
    if (x->currentFrame>0) 
    {
        for(i=0; i<x->featureLength; i++)
        {
            sum = 0;
            for(j=0; j<x->currentFrame; j++) 
            {
                sum = sum + x->instances[j].instance[i];
            }
            sum = sum / x->currentFrame;
            atom_setfloat(x->x_listOut+i, sum); 
        }
        
        outlet_list(x->featureList, 0, x->featureLength, x->x_listOut);
        
        x->currentFrame = 0;
    
    }
}

static void featureMean_clear(t_featureMean *x)
{
	int i, j;
    
	// free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
	t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
    
	x->currentFrame = 0;
    
    for(i=0; i<x->featureLength; i++)
        atom_setfloat(x->x_listOut+i, 0.0);
    
    x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
    
	for(i=0; i<x->numFrames; i++)
		x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
    
	for(i=0; i<x->numFrames; i++)
		for(j=0; j<x->featureLength; j++)
			x->instances[i].instance[j] = 0.0;
}

static void featureMean_numFrames(t_featureMean *x, int num)
{
	int i, j;
    
	if(num)
	{
        //x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, x->featureLength*sizeof(t_atom), (x->featureLength*num)*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->currentFrame = 0;
		x->numFrames = num;
        
        for(i=0; i<x->featureLength; i++)
	        atom_setfloat(x->x_listOut+i, 0.0);
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureMean_length(t_featureMean *x, int len)
{
	int i, j;
    
	if(len)
	{
        x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, x->featureLength*sizeof(t_atom), len*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		x->featureLength = len;
		x->currentFrame = 0;
        
        for(i=0; i<x->featureLength; i++)
	        atom_setfloat(x->x_listOut+i, 0.0);
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureMean_print(t_featureMean *x)
{
    post("averaging %i vectors with %i features", x->numFrames, x->featureLength);
}

static void featureMean_free(t_featureMean *x)
{
	int i;      
    
    // free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
    t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
    
    // free listOut memory
	t_freebytes_(x->x_listOut, x->featureLength*sizeof(t_atom));
    
}

static void *featureMean_new(t_symbol *s, long argc, t_atom *argv)
{
    t_featureMean *x = NULL;
    int i,j;
    t_atom *ap;

	if ((x = (t_featureMean *)object_alloc(featureMean_class))) {
        
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
     
        x->x_listOut = (t_atom *)t_getbytes_(x->featureLength*sizeof(t_atom));
        x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
        for(i=0; i<x->featureLength; i++)
            atom_setfloat(x->x_listOut+i, 0.0);
        
        for(i=0; i<x->numFrames; i++)
            x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
        for(i=0; i<x->numFrames; i++)
            for(j=0; j<x->featureLength; j++)
                x->instances[i].instance[j] = 0.0;

    }
    
    return (x);
}

