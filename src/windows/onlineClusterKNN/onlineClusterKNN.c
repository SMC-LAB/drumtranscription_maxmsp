#include "ext.h"							// standard Max include, always required
#include "ext_obex.h"						// required for new style Max object
#include "m_memory.h"	
#define MAXCLUSTERMEMS 8192
#define t_float float

////////////////////////// object struct
typedef struct instance
{
    float *instance;
    int cluster;
} t_instance;

typedef struct _onlineClusterKNN
{
	t_object  x_obj;
    t_instance *means; //means for each cluster
    int *num; //counter of instances for each cluster
    
    t_instance *instances;
	t_float *featureInput;
    int *instanceFeatureLengths;
    int featureLength;
    int numInstances;
    int numClusters;
    int distMetric; 
    int randomFlag;
    int memory;
    
    void *id;
} t_onlineClusterKNN;


///////////////////////// function prototypes
//static void printv(t_float *v, int listLength);
int countDigits(int num);
static int test_zero(t_float *newelement, int numfeatures);
static void compute_mean(t_onlineClusterKNN *x, int k, t_instance newelement);
static t_float squared_euclid(t_float *v1, t_float *v2, int numfeatures);

static void onlineClusterKNN_memory(t_onlineClusterKNN *x, int mem);
static void onlineClusterKNN_print(t_onlineClusterKNN *x);
static void onlineClusterKNN_clear(t_onlineClusterKNN *x);
static void onlineClusterKNN_cluster(t_onlineClusterKNN *x, t_symbol *s, int argc, t_atom *argv);
//static void onlineClusterKNN_bang(t_onlineClusterKNN  *x);
static void onlineClusterKNN_random(t_onlineClusterKNN  *x, int r);
static void onlineClusterKNN_free(t_onlineClusterKNN *x);
static void *onlineClusterKNN_new(t_symbol *s, long argc, t_atom *argv);

//////////////////////// global class pointer variable
void *onlineClusterKNN_class;


int main(void)
{	
	t_class *c;
	
	c = class_new("onlineClusterKNN", (method)onlineClusterKNN_new, (method)onlineClusterKNN_free, (long)sizeof(t_onlineClusterKNN), 
				  0L /* leave NULL!! */,  A_GIMME, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	onlineClusterKNN_class = c;
    
    
    class_addmethod(onlineClusterKNN_class, (method)onlineClusterKNN_cluster, "list", A_GIMME, 0);   
    
    class_addmethod(onlineClusterKNN_class, 
                    (method)onlineClusterKNN_clear,
                    "clear",
                    0
                    );  
    
    class_addmethod(onlineClusterKNN_class, 
                    (method)onlineClusterKNN_print,
                    "print",
                    0
                    );
    
    class_addmethod(onlineClusterKNN_class, 
                    (method)onlineClusterKNN_random,
                    "random",
                    A_LONG,
                    0
                    );
    
    class_addmethod(onlineClusterKNN_class, 
                    (method)onlineClusterKNN_memory,
                    "memory",
                    A_LONG,
                    0
                    );

	post("onlineClusterKNN version 0.1.0");
	return 0;
}


static void onlineClusterKNN_print(t_onlineClusterKNN *x)
{
    //int i,j;
	post("no. of instances: %i", x->numInstances);
	post("feature length: %i", x->featureLength);	
	post("distance metric: %i", x->distMetric);
	post("no. of clusters: %i", x->numClusters);
    post("memory: %i", x->memory);
    post("initialize with random means: %i", x->randomFlag);
	/*
    if (x->numInstances > x->numClusters - 1)
        for (i=0; i<x->numClusters; i++) {
            post("mean %i = ", i);
            for(j=0; j<x->featureLength; j++)
            {
                post("%f ",x->means[i].instance[j]);
            }
        }
     */
}

static void onlineClusterKNN_clear(t_onlineClusterKNN *x)
{
	int i,j;	
    
    if (x->randomFlag == 1)
    {
        // free the database memory
        for(i=0; i<x->numInstances; i++)
            t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
        
        for(i=0; i<x->numClusters; i++)
            t_freebytes_(x->means[i].instance, x->instanceFeatureLengths[0]*sizeof(t_float));
        
        x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), 0);
        x->means = (t_instance *)t_resizebytes_(x->means, x->numClusters * sizeof(t_instance), 0);
        x->num = (int *)t_resizebytes_(x->num, x->numClusters * sizeof(int), 0);
        x->numInstances = 0;
        x->featureLength = 0;
        x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), 0);	
    }
    else
    {              
        for(i=0; i<x->numClusters; i++)
        {
            for(j=0; j<x->featureLength; j++)
            {
                x->means[i].instance[j] =  x->instances[i].instance[j];
            }
            x->num[i] = 1;
        }
        
        // free the database memory
        for(i=x->numClusters; i<x->numInstances; i++)
            t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
        
        x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), x->numClusters * sizeof(t_instance));
        x->numInstances = x->numClusters;
        x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), x->numClusters * sizeof(int));	        
    }   
    
    post("all instances cleared.");
}


static void onlineClusterKNN_cluster(t_onlineClusterKNN *x, t_symbol *s, int argc, t_atom *argv)
{
    int i, j, k, instanceIdx, listLength;
    float min_dist,dist;
    
	instanceIdx = x->numInstances;
	listLength = argc;
	//s=s; // to get rid of 'unused variable' warning
    //post("list length: %i", listLength);
    
    if((x->featureLength>0) && (x->featureLength != listLength))
	{
        post("received list of length %i and expected %i", listLength, x->featureLength); 
        return;
    }    
    
	x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), (x->numInstances+1) * sizeof(t_instance));
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), (x->numInstances+1) * sizeof(int));
	
	x->instanceFeatureLengths[instanceIdx] = listLength;
	x->instances[instanceIdx].instance = (t_float *)t_getbytes_(listLength * sizeof(t_float));    
    
	x->numInstances++;
	//post("no. of instances: %i", x->numInstances);          
    
    
    //get the data
	for(i=0; i<listLength; i++)
		x->instances[instanceIdx].instance[i] = atom_getfloat(argv+i);
    
    //test if received element is zeros vector
    if (test_zero(x->instances[instanceIdx].instance,listLength) == 0)
    {
        //post("instance cannot be zeros vector");
        //rollback
        t_freebytes_(x->instances[instanceIdx].instance, x->instanceFeatureLengths[instanceIdx]*sizeof(t_float));
        x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), (x->numInstances-1) * sizeof(t_instance));
        x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), (x->numInstances-1) * sizeof(int)); 
        
        x->numInstances--;
        
    }
    else 
    {
        
        if(instanceIdx == 0)
        {
            x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));		
            x->featureLength = listLength;  
            
            x->num = (int *)t_resizebytes_(x->num, 0 * sizeof(int), x->numClusters * sizeof(int));
            
            // initialize means randomly for each cluster        
            for(i=0; i<x->numClusters; i++)
            {
                x->means = (t_instance *)t_resizebytes_(x->means, i * sizeof(t_instance), (i+1) * sizeof(t_instance));
                x->means[i].instance = (t_float *)t_getbytes_(listLength * sizeof(t_float)); 
                
                if (x->randomFlag == 1)
                {
                    for(j=0; j<listLength; j++)
                    {
                        srand(i*j+i+j+1);
                        x->means[i].instance[j] = (float)rand()/(float)RAND_MAX;            
                    }
                }
            }                
            // initialize number of instances for each cluster
            for(i=0; i<x->numClusters; i++)
                x->num[i] = 0;
        };
        
        
        //normalize the data to be 0-1
        for(i=0; i<listLength; i++)
            if (x->instances[instanceIdx].instance[i]>1)
            {            
                x->instances[instanceIdx].instance[i] = x->instances[instanceIdx].instance[i] * pow(10,-countDigits((int)(x->instances[instanceIdx].instance[i])));
            }
        
        //////////////ONLINE CLUSTERING
        
        //initialize means with the first instances if random==0
        if ((x->randomFlag == 0) && (instanceIdx < x->numClusters))
        {
            for(j=0; j<listLength; j++)
            {
                x->means[instanceIdx].instance[j] =  x->instances[instanceIdx].instance[j];
            }
            x->instances[instanceIdx].cluster = instanceIdx;
            x->num[instanceIdx] = 1;
        }
        else
        {        
            //compute distances to the means and determine the closest cluster 
            min_dist = 99999999;
            k = -1;
            for(i=0; i<x->numClusters; i++)
            { 
                dist = squared_euclid(x->means[i].instance,x->instances[instanceIdx].instance,listLength);
                //post("%6.20f", dist);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    k = i;
                }
            }   
            
            //add the new instance to the found cluster
            //post("cluster %i", k);
            if (k != -1)
            {
                x->num[k] = x->num[k] + 1;
                compute_mean(x, k, x->instances[instanceIdx]);
                x->instances[instanceIdx].cluster = k;     
            }  
        }
        
        outlet_int(x->id, x->instances[instanceIdx].cluster);    
    }
}


static void onlineClusterKNN_random(t_onlineClusterKNN  *x, int r)
{
	r = (r<0)?0:r;
	r = (r>1)?1:r;
	x->randomFlag = r;
    post("randomize initial means: %i", x->randomFlag);
}

static void onlineClusterKNN_memory(t_onlineClusterKNN  *x, int m)
{
	m = (m<0)?0:m;
	x->memory = m;
    post("online memory: %i", x->memory);
}

static void onlineClusterKNN_free(t_onlineClusterKNN *x)
{
	int i;
    i = 0;        
    
	// free the database memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
    
    for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->means[i].instance, x->instanceFeatureLengths[0]*sizeof(t_float));
    
    t_freebytes_(x->means, x->numClusters*sizeof(t_instance));
    t_freebytes_(x->num, x->numClusters*sizeof(int));    
	
	t_freebytes_(x->instances, x->numInstances*sizeof(t_instance)); 
	
	t_freebytes_(x->featureInput, x->featureLength*sizeof(t_float));
	t_freebytes_(x->instanceFeatureLengths, x->numInstances*sizeof(int));	
    
}

static void *onlineClusterKNN_new(t_symbol *s, long argc, t_atom *argv)
{
    t_onlineClusterKNN *x = NULL;
    
    // object instantiation, NEW STYLE
	if ((x = (t_onlineClusterKNN *)object_alloc(onlineClusterKNN_class))) {

        
        if (argc > 0 && (argv)->a_type == A_LONG){
            //is first arg an int - then set to numClusters
            x->numClusters = atom_getlong(argv); 
        } else {
            x->numClusters = 2;
        }
        object_post((t_object *)x, "number of clusters %ld", x->numClusters);
        
        x->id = intout(x);
        //x->id = outlet_new (x, NULL);
        
        x->means = (t_instance *)t_getbytes_(0);
        x->num = (int *)t_getbytes_(0);
        
        x->instances = (t_instance *)t_getbytes_(0);
        x->instanceFeatureLengths = (int *)t_getbytes_(0);
        x->featureInput = (t_float *)t_getbytes_(0);
        
        x->featureLength = 0;
        x->memory = 0;
        x->numInstances = 0;
        x->distMetric = 0;  // euclidean distance by default
        x->randomFlag = 0;  //do not initialize the knn means with random values by default
        x->numClusters = (x->numClusters>10)?10:x->numClusters;
        x->numClusters = (x->numClusters<2)?2:x->numClusters;
    }
    
    return (x);
}


/* ---------------- utility functions ---------------------- */

int countDigits(int num){
    int count=0;
    
    while(num>0){
        num=num/10;
        count++;
    }
    
    return count;
}

static int test_zero(t_float *newelement, int numfeatures)
{
	int i,b;
    b=0;
	for(i=0; i<numfeatures; i++)
	{		
        if (newelement[i]!=0) {b=1;break;}
    }            
    return(b);
}

static void compute_mean(t_onlineClusterKNN *x, int k, t_instance newelement)
{
	int i;
    
	for(i=0; i<x->featureLength; i++)
	{	
        if (x->memory == 0)
        x->means[k].instance[i] = x->means[k].instance[i] + (newelement.instance[i] - x->means[k].instance[i])/x->num[k];
        else x->means[k].instance[i] = x->means[k].instance[i] + (newelement.instance[i] - x->means[k].instance[i])/x->memory;
	}
}

static t_float squared_euclid(t_float *v1, t_float *v2, int numfeatures)
{
	int i;
	t_float sum, dist;
	
	sum=dist=0;
	for(i=0; i < numfeatures; i++)
	{
        sum = sum + pow((v1[i]-v2[i]),2.0);	
    }
    dist = sqrt(sum);
	return(dist);
}

/* ---------------- END utility functions ---------------------- */