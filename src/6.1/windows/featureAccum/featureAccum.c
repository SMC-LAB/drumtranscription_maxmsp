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

typedef struct _featureAccum
{
	t_object  x_obj;
    t_atom *x_listOut;
    t_instance *instances;
    int spew;
    int featureLength;
    int currentFrame;
    int numInstances;
    int numFrames;
    
    void *featureList;
} t_featureAccum;


///////////////////////// function prototypes
static void featureAccum_accum(t_featureAccum *x, t_symbol *s, int argc, t_atom *argv);
static void featureAccum_numFrames(t_featureAccum *x, int num);
static void featureAccum_spew(t_featureAccum *x, int s);
static void featureAccum_length(t_featureAccum *x,int len);
static void featureAccum_clear(t_featureAccum *x);
static void featureAccum_free(t_featureAccum *x);
static void *featureAccum_new(t_symbol *s, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *featureAccum_class;


int main(void)
{	
	t_class *c;
	
	c = class_new("featureAccum", (method)featureAccum_new, (method)featureAccum_free, (long)sizeof(t_featureAccum),0L /* leave NULL!! */,  A_GIMME, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	featureAccum_class = c;
        
    class_addmethod(featureAccum_class, (method)featureAccum_accum, "list", A_GIMME, 0);
    
    
    
    class_addmethod(featureAccum_class, 
                    (method)featureAccum_clear,
                    "clear",
                    0
                    );
    
    class_addmethod(featureAccum_class, 
                    (method)featureAccum_numFrames,
                    "num_frames",
                    A_LONG,
                    0
                    );
    
    class_addmethod(featureAccum_class, 
                    (method)featureAccum_spew,
                    "spew",
                    A_LONG,
                    0
                    );
    
    class_addmethod(featureAccum_class, 
                    (method)featureAccum_length,
                    "length",
                    A_LONG,
                    0
                    );


	post("featureAccum version 0.1.0");
	return 0;
}

static void featureAccum_accum(t_featureAccum *x, t_symbol *s, int argc, t_atom *argv)
{
	int i, j, count, totalFeat;
    
	if(x->featureLength != argc)
		error("featureAccum: input length does not match current length setting. input ignored.");
	else
		for(i=0; i<x->featureLength; i++)
			x->instances[x->currentFrame].instance[i] = atom_getfloat(argv+i);
    
	x->currentFrame++;
    
	if((x->currentFrame==x->numFrames) || (x->spew))
	{
		totalFeat = x->featureLength * x->numFrames;
        
		for(count=0, i=x->numFrames-x->currentFrame; count<x->numFrames; count++, i++)
			for(j=0; j<x->featureLength; j++)
				atom_setfloat(x->x_listOut+((i%x->numFrames)*x->featureLength)+j, x->instances[count].instance[j]);
        
		outlet_list(x->featureList, 0, totalFeat, x->x_listOut);
        
		x->currentFrame = (x->currentFrame==x->numFrames)?0:x->currentFrame;
	}
}

static void featureAccum_clear(t_featureAccum *x)
{
	int i, j;
    
	// free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
	t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
    
	x->currentFrame = 0;
    
    for(i=0; i<x->featureLength*x->numFrames; i++)
        atom_setfloat(x->x_listOut+i, 0.0);
    
    x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
    
	for(i=0; i<x->numFrames; i++)
		x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
    
	for(i=0; i<x->numFrames; i++)
		for(j=0; j<x->featureLength; j++)
			x->instances[i].instance[j] = 0.0;
}

static void featureAccum_numFrames(t_featureAccum *x, int num)
{
	int i, j;
    
	if(num)
	{
        x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, (x->featureLength*x->numFrames)*sizeof(t_atom), (x->featureLength*num)*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->currentFrame = 0;
		x->numFrames = num;
        
        for(i=0; i<x->featureLength*x->numFrames; i++)
	        atom_setfloat(x->x_listOut+i, 0.0);
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureAccum_length(t_featureAccum *x, int len)
{
	int i, j;
    
	if(len)
	{
        x->x_listOut = (t_atom *)t_resizebytes_(x->x_listOut, (x->featureLength*x->numFrames)*sizeof(t_atom), (len*x->numFrames)*sizeof(t_atom));
        
		// free the database memory
		for(i=0; i<x->numFrames; i++)
			t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
        
		t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
		x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
		x->featureLength = len;
		x->currentFrame = 0;
        
        for(i=0; i<x->featureLength*x->numFrames; i++)
	        SETFLOAT(x->x_listOut+i, 0.0);
        
		for(i=0; i<x->numFrames; i++)
			x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
		for(i=0; i<x->numFrames; i++)
			for(j=0; j<x->featureLength; j++)
				x->instances[i].instance[j] = 0.0;
	}
}

static void featureAccum_spew(t_featureAccum *x, int s)
{
	s = (s<=0)?0:s;
	s = (s>=1)?1:s;
	x->spew = s;
    
	post("spew mode: %i", x->spew);
}



static void featureAccum_free(t_featureAccum *x)
{
	int i;      
    
    // free the database memory
	for(i=0; i<x->numFrames; i++)
		t_freebytes_(x->instances[i].instance, x->featureLength*sizeof(float));
    
    t_freebytes_(x->instances, x->numFrames*sizeof(t_instance));
        
    // free listOut memory
	t_freebytes_(x->x_listOut, x->featureLength*sizeof(t_atom));
    
}

static void *featureAccum_new(t_symbol *s, long argc, t_atom *argv)
{
    t_featureAccum *x = NULL;
    int i,j;
    t_atom *ap;

	if ((x = (t_featureAccum *)object_alloc(featureAccum_class))) {
        
        x->featureLength = 2;
        x->numFrames = 2;
        x->currentFrame = 0;
        x->spew = 0;
        // increment ap each time to get to the next atom
        for (i = 0, ap = argv; i < argc; i++, ap++) 
        {
            if (atom_gettype(ap) == A_LONG) 
            {
                if (i==0) x->featureLength = atom_getlong(ap);
                if (i==1) x->numFrames = atom_getlong(ap);
                if (i==2) x->spew = (atom_getlong(ap) < 0) ? 0 : (atom_getlong(ap) > 1) ? 1 : atom_getlong(ap);
            }
            else post("%ld: unknown atom type (%ld)", i+1, atom_gettype(ap));
        }
        
        
        x->featureList = listout(x);
        
        x->x_listOut = (t_atom *)t_getbytes_((x->featureLength*x->numFrames)*sizeof(t_atom));
        x->instances = (t_instance *)t_getbytes_(x->numFrames*sizeof(t_instance));
        
        for(i=0; i<x->featureLength*x->numFrames; i++)
            atom_setfloat(x->x_listOut+i, 0.0);
        
        for(i=0; i<x->numFrames; i++)
            x->instances[i].instance = (float *)t_getbytes_(x->featureLength*sizeof(float));
        
        for(i=0; i<x->numFrames; i++)
            for(j=0; j<x->featureLength; j++)
                x->instances[i].instance[j] = 0.0;

    }
    
    return (x);
}

