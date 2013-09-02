/*
  THIS CODE WAS ADAPTED FOR MAXMSP FROM THE TIMBREID PURE DATA LIBRARY BY MARIUS MIRON, 2012, SMC GROUP, INESC PORTO, PORTUGAL
 timbreID - A generic classification external.
 
 Copyright 2009 William Brent
 
 This file is part of timbreID.
 
 timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
 timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 
 version 0.6.0, December 23, 2011
 
 • 0.6.0 changed all floats to t_floats, got rid of unused headers, removed underscores from variable names.  FIXED concatenative_id function's search center/neighborhood bug.  Added a similarity matrix function.  Added a "print" method to show internal settings.  Made even or odd neighborhood settings valid since it will sometimes be important to pick an exact range.  Odd settings are fine - with n=5, we search from the center -/+ 2.  With n=6, it be center -2/+3 because of a special case for searchStart.
 • 0.5.1 was just to remove the "no such table" error from the NRT externs when it doesn't find the specified table at load time. also, the max window size for all the NRT externs was set to 131072 because this seems to be the maximum allowable size for mayer_realfft().
 • 0.5 as part of the timbreID-0.5 update, this gets rid of unnecessary getbytes(0) calls in timbreID.c. does not yet use memset() or memcpy() where possible.  fixed lack of memory free for x->weights, x->attributeOrder, x->instanceFeatureLengths, x->normData in _free.
 • fixed small memory mismanagement in all _make_filterbank functions for bark-based externs.
 • [timbreID] is unchanged - only updating for addition of [featureAccum] and [binWrangler], which streamline time-evolving feature extraction
 • [timbreID] is unchanged - only updating version for [cepstrum] and [cepstrum~] changes.
 • part of the update acknowledging the change to one tID folder full of sources and binaries. This update also included filter power averaging in barkSpec, bfcc, and mfcc.
 • removed cosine similarity distance metric
 • 0.3.3 reflects the general package update, which includes specIrregularity(~), and fixes a bug in cepstrum(~). also, changed timbreID's default KNN setting to 1.  It should really only be anything different if clustering is done.
 • 0.3.3 is based on the windows source, where the _read functions were changed to avoid using fdopen.  also corrected the use of a resizable array for semicolon in the _read_cluster_text function. also added static to all functions except _setup.
 • 0.3.2 no changes to [timbreID]. Added RT and NRT spectral spread externs, plus mean and standard deviation externs for summarizing spectro-temporal features of various lengths
 • 0.3.1G this gets around to adding the cluster read/write functions. F is cleaned up version of D. 0.3.1G does work with read/write clusters.
 • 0.3.1F skips the different approach in E (didn't take that very far), and fixed D to work.  The problem in _compute_cluster was after clustering and in memory cleanup. Shifting elements of x->clusterMembers and x->clusterMemberLengths was done in the wrong order and memory was resized incorrectly. F cleans up D, and adds the cluster_write functions.
 • 0.3.1D all that remains is cluster write/read.
 • 0.3.1C extends 0.3.1B to shoot for making x->clusterMembers memory dynamically allocated as well.  C works up to and including unclustering - next will be the read/write functions.
 • 0.3.1 finally makes feature database and cluster member memory dynamically allocated. Set up the _read and _write functions so that the header includes the entire instanceFeatureLengths array. Might as well keep all that information so that mixed instance databases won't be a problem in the future. BUT - write/read_text do not store/retrieve this information. More importantly, x->featureLength is simply set to the first element in x->instanceFeatureLengths in the _read function. So, the repercussions of mixed length features in the database has not been thought through. A common feature length is still assumed.
 • 0.3.0 adds a method for finding the worst match, and two others for outputting the maxes and mins of all features in the database (as long as it has been normalized). Also - found a bug in the concat_id chain. timbreID_clear never set x->neighborhood back to 0, which caused crashes after training/clearing/training a few times.
 • 0.2.9 removes some warnings about unitialized variables. the entire timbreID set of externs was updated on this date as well to include non real-time versions of all feature externs. 0.3.0 will finally implement dynamic memory allocation for feature database instances, rather than using MAXFEATURELENGTH for all (very wasteful).
 • 0.2.8 feature externs now use mayer_realfft, switched to timbreID-help.pd rather than help-timbreID.pd for the help patch name.
 • 0.2.7C Added clusterMembership function
 • 0.2.7 Added ability to output a instance's feature list from the 4th outlet.
 • 0.2.6 Added ability to write and read cluster instances in binary and text formats.
 • 0.2.5 Added correlation as a distance metric. Fixed a bug in compute_cluster: it now throws an error if you ask for 0 clusters.  Also slightly tweaked the clustering algorithm itself.
 • 0.2.4 Fixed major bug where timbreID_train wasn't updating the size of x->featureInput if the feature length changes. This caused crashes in timbreID_id after training with a feature length different than the default of 47 points. Fixed bug where requesting compute_cluster with a different number of clusters when already clustered causes crash. Fixed error in _write: j was used uninitialized. Distance metric choice is now used in compute_cluster. Previously, it was squared_euclid only. Added method for outputting cluster member lists, or the entire cluster list. 
 • 0.2.3 adds binary file output (.timid now default), and an option for text output.
 • 0.2.2 adds MATLAB .mat file output.
 • 0.2.0 adds cosine distance function (doesn't account for attribute normalization yet though!).
 • 0.1.9 adds ARFF function to export data for WEKA.
 
 • fixed bug where the distances of the final elements in a compute_order result were sometimes 
 larger than INT_MAX, so that they were all index 0 (init value). Using FLT_MAX instead now, which 
 should be big enough to be safe.
 
 • Accounted for normalization flag in the compute_variance function (requiring a new timbreID_mean 
 function as well).
 
 */


#include "ext.h"							// standard Max include, always required
#include "ext_obex.h"						// required for new style Max object
#include "m_memory.h"
#include <math.h>
#include <float.h>
#define MAXCLUSTERMEMS 8192
#define MAXPDSTRING 1000        /* must be >= FILENAME_MAX */
#define MAXPDARG 5              /* max number of args we can typecheck today */
#define t_float float

////////////////////////// object struct

typedef struct instance
{
    t_float *instance;
} t_instance;

typedef struct member
{
    int *member;
} t_member;

typedef struct knn_info
{
    t_float dist;
	t_float safeDist;
    int idx;
    int cluster;
} t_knn_info;

typedef struct normData
{
    t_float max;
    t_float min;
	t_float denominator;
} t_normData;


typedef struct _timbreID
{
    t_object x_obj;
    t_instance *instances;
    t_member *clusterMembers;
    int *clusterMemberLengths;
    t_knn_info *knnDistsIdxs;
    t_normData *normData;
	t_float *featureInput;
    int *instanceClusterMembership;
    int *instanceFeatureLengths;
    int featureLength;
    int numInstances;
    int numClusters;
    int distMetric;
    int k;
    int normalize;
    int relativeOrdering;
    // must deal with resizing attributeOrder and weights next...
    int *attributeOrder;
    t_float *weights;
    
    int reorientFlag;
    int neighborhood;
    int searchCenter;
    int prevMatch;
    int maxMatches;
    int stutterProtect;
    t_float jumpProb;
    
    int attributelo;
    int attributehi;
    
    void *id;
    void *nearestDist;
    void *confidence;
    void *x_orderList;
    void *x_featureList;
    
    void		*proxy1, *proxy2;			// proxy inlet
    long		proxy_inletnum1, proxy_inletnum2;	// # of inlet currently in use
    
    //t_canvas *x_canvas;
   
} t_timbreID;



///////////////////////// function prototypes
static void timbreID_makepath(char *filename, char *buf);
static void timbreID_sort_knn_info(int k, int numInstances, int prevMatch, t_knn_info *list);
static void timbreID_sort_float(int n, t_float *list);
static t_float timbreID_mean(int num_rows, int column, t_instance *instances, int normal_flag, t_normData *normData);
static t_float timbreID_squared_euclid(t_timbreID *x, t_float *v1, t_float *v2);
static t_float timbreID_manhattan(t_timbreID *x, t_float *v1, t_float *v2);
static t_float timbreID_correlation(t_timbreID *x, t_float *v1, t_float *v2);

static void timbreID_input(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_train(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_worst_match(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_forget(t_timbreID *x);

static void timbreID_concat_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_concat_neighborhood(t_timbreID *x, int n);
static void timbreID_concat_jumpProb(t_timbreID *x,int jp);
static void timbreID_concat_reorient(t_timbreID *x, int r);
static void timbreID_concat_searchCenter(t_timbreID *x, int sc);
static void timbreID_concat_maxMatches(t_timbreID *x, int mm);
static void timbreID_concat_stutterProtect(t_timbreID *x, int sp);

static void timbreID_knn(t_timbreID *x, int k);
static void timbreID_normalize(t_timbreID *x, int n);
static void timbreID_print(t_timbreID *x);
static void timbreID_numInstances(t_timbreID *x);
static void timbreID_manual_cluster(t_timbreID *x, int numClusters, int cluster_idx, int low, int hi);
static void timbreID_compute_cluster(t_timbreID *x, int numClusters);
static void timbreID_uncluster(t_timbreID *x);
static void timbreID_compute_variance(t_timbreID *x);
static void timbreID_clusters_list(t_timbreID *x);
static void timbreID_cluster_list(t_timbreID *x, int idx);
static void timbreID_clusterMembership(t_timbreID *x, int idx);
static void timbreID_compute_order(t_timbreID *x, int reference);
static void timbreID_relativeOrdering(t_timbreID *x, int rel);
static void timbreID_distMetric(t_timbreID *x, int f);
static void timbreID_weights(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_attributes(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_attribute_range(t_timbreID *x, int lo, int hi);
static void timbreID_order_attributes(t_timbreID *x);
static void timbreID_print_instance(t_timbreID *x, int e, int f, int g);
static void timbreID_feature_list(t_timbreID *x, int idx);
static void timbreID_similarityMatrix(t_timbreID *x, int startInstance, int finishInstance, int normalize);
static void timbreID_max_values(t_timbreID *x);
static void timbreID_min_values(t_timbreID *x);

static void timbreID_clear(t_timbreID *x);
static void timbreID_write(t_timbreID *x, t_symbol *s);
static void timbreID_read(t_timbreID *x, t_symbol *s);
static void timbreID_write_text(t_timbreID *x, t_symbol *s);
static void timbreID_read_text(t_timbreID *x, t_symbol *s);
static void timbreID_ARFF(t_timbreID *x, t_symbol *s, int argc, t_atom *argv);
static void timbreID_MATLAB(t_timbreID *x, t_symbol *file_symbol, t_symbol *var_symbol);
static void timbreID_write_clusters(t_timbreID *x, t_symbol *s);
static void timbreID_read_clusters(t_timbreID *x, t_symbol *s);
static void timbreID_write_clusters_text(t_timbreID *x, t_symbol *s);
static void timbreID_read_clusters_text(t_timbreID *x, t_symbol *s);

static void timbreID_free(t_timbreID *x);
static void *timbreID_new(t_symbol *s);



//////////////////////// global class pointer variable
void *timbreID_class;


int C74_EXPORT main(void)
{	
	t_class *c;
	
	c = class_new("timbreID", (method)timbreID_new, (method)timbreID_free, (long)sizeof(t_timbreID), 
				  0L /* leave NULL!! */,  A_DEFFLOAT, 0);

	
	class_register(CLASS_BOX, c); /* CLASS_NOBOX */
	timbreID_class = c;
    
    //class_addmethod(timbreID_class, (method)timbreID_train, "list", A_GIMME, 0);   
    //class_addmethod(timbreID_class, (method)timbreID_id, "in1", A_GIMME, 0);
    //class_addmethod(timbreID_class, (method)timbreID_concat_id, "in2", A_GIMME, 0);
    
    class_addmethod(timbreID_class, (method)timbreID_input, "list", A_GIMME, 0); 
    

	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_worst_match,
                    "worst_match",
                    A_GIMME,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_forget, 
                    "forget",
                    0
                    );

	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_neighborhood,
                    "neighborhood",
                     A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_jumpProb,
                    "jump_prob",
                     A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_reorient,
                    "reorient",
                     A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_searchCenter,
                    "search_center",
                     A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_maxMatches,
                    "max_matches",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_concat_stutterProtect,
                    "stutter_protect",
                    A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_knn, 
                    "knn",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_normalize, 
                    "normalize",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_print,
                    "print",
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_numInstances,
                    "num_instances",
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_manual_cluster, 
                    "manual_cluster",
                    A_DEFLONG,
                    A_DEFLONG,
                    A_DEFLONG,
                    A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_compute_cluster, 
                    "cluster",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_uncluster, 
                    "uncluster",
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_clusters_list, 
                    "clusters_list",
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_cluster_list, 
                    "cluster_list",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_clusterMembership, 
                    "cluster_membership",
                    A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_compute_order, 
                    "order",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_relativeOrdering, 
                    "relative_ordering",
                    A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_compute_variance, 
                    "variance",
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_distMetric,
                    "dist_metric",
                    A_DEFLONG,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_weights,
                    "weights",
                    A_GIMME,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_attributes,
                    "attributes",
                    A_GIMME,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_attribute_range,
                    "attribute_range",
                    A_DEFLONG,
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_order_attributes, 
                    "order_attributes",
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_print_instance,
                    "print_instance",
                    A_DEFLONG,
                    A_DEFLONG,
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_feature_list, 
                    "feature_list",
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_similarityMatrix, 
                    "similarity_matrix",
                    A_DEFLONG,
                    A_DEFLONG,
                    A_DEFLONG,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_max_values, 
                    "max_values",
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_min_values, 
                    "min_values",
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_clear,
                    "clear",
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_write,
                    "write",
                    A_SYM,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_read,
                    "read",
                    A_SYM,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_write_text,
                    "write_text",
                    A_SYM,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_read_text,
                    "read_text",
                    A_SYM,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_ARFF,
                    "ARFF",
                    A_GIMME,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_MATLAB,
                    "MATLAB",
                    A_SYM,
                    A_SYM,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_write_clusters,
                    "write_clusters",
                    A_SYM,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_read_clusters,
                    "read_clusters",
                    A_SYM,
                    0
                    );
	
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_write_clusters_text,
                    "write_clusters_text",
                    A_SYM,
                    0
                    );
    
	class_addmethod(
                    timbreID_class, 
                    (method)timbreID_read_clusters_text,
                    "read_clusters_text",
                    A_SYM,
                    0
                    );


	post("timbreID version 0.6.0");
	return 0;
}


static void timbreID_input(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
    long inlet = proxy_getinlet((t_object *)x); // what inlet did this message come in through?

    switch (inlet) {
        case 1:
            //post("right");
            timbreID_concat_id(x,s,argc,argv);
            break;
        case 2:
            //post("middle");
            timbreID_id(x,s,argc,argv);
            break;    
        default:
            //post("left");
            timbreID_train(x,s,argc,argv);
            break;
    }

}


static void *timbreID_new(t_symbol *s)
{
    t_timbreID *x = NULL;
    
    // object instantiation, NEW STYLE
	if ((x = (t_timbreID *)object_alloc(timbreID_class))) {
        //intin(x, 1);
        //intin(x, 2);
        x->proxy1 = proxy_new(x, 1, &x->proxy_inletnum1);	// fully-flexible inlet for any type 
        x->proxy2 = proxy_new(x, 2, &x->proxy_inletnum2);	// fully-flexible inlet for any type         
        
        x->x_featureList = listout((t_object *)x); 
        x->x_orderList = listout((t_object *)x);
        x->confidence = floatout((t_object *)x);
        x->nearestDist = floatout((t_object *)x);
        x->id = intout((t_object *)x);        
        
        
        x->instances = (t_instance *)t_getbytes_(0);
        x->instanceFeatureLengths = (int *)t_getbytes_(0);
        x->clusterMembers = (t_member *)t_getbytes_(0);
        x->clusterMemberLengths = (int *)t_getbytes_(0);
        x->instanceClusterMembership = (int *)t_getbytes_(0);
        x->weights = (t_float *)t_getbytes_(0);
        x->attributeOrder = (int *)t_getbytes_(0);
        x->knnDistsIdxs = (t_knn_info *)t_getbytes_(0);
        x->featureInput = (t_float *)t_getbytes_(0);
        
        x->featureLength = 0;
        x->numClusters=0;
        x->numInstances = 0;
        x->distMetric = 0;  // euclidean distance by default
        x->k = 1;
        x->normalize = 0;
        x->relativeOrdering = 1;
        x->stutterProtect = 0;
        
        x->prevMatch = -1;
        x->maxMatches = 3;
        x->reorientFlag = 0;
        x->neighborhood = 0;
        x->searchCenter = 0;
        x->jumpProb = 0.0;
    }
    
    return (x);
}

static void timbreID_free(t_timbreID *x)
{
	int i;
    
	// free the database memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
	
	t_freebytes_(x->instances, x->numInstances*sizeof(t_instance));
    
	// free the clusterMembers memory
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
	t_freebytes_(x->clusterMembers, x->numClusters*sizeof(t_member));
    
	// free cluster member lengths
	t_freebytes_(x->clusterMemberLengths, x->numClusters*sizeof(int));
    
	t_freebytes_(x->knnDistsIdxs, x->numInstances*sizeof(t_knn_info));
	t_freebytes_(x->instanceClusterMembership, x->numInstances*sizeof(int));
	t_freebytes_(x->featureInput, x->featureLength*sizeof(t_float));
	t_freebytes_(x->weights, x->featureLength*sizeof(t_float));
	t_freebytes_(x->attributeOrder, x->featureLength*sizeof(int));
	t_freebytes_(x->instanceFeatureLengths, x->numInstances*sizeof(int));
	if(x->normalize)
		t_freebytes_(x->normData, x->featureLength*sizeof(t_normData));
    
    object_free(x->proxy1);		// frees all resources associated with the proxy
    object_free(x->proxy2);	
}

/* ------------------------ timbreID -------------------------------- */

static void timbreID_train(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i, instanceIdx, listLength;
    
	instanceIdx = x->numInstances;
	listLength = argc;
	//s=s; // to get rid of 'unused variable' warning
    
	x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), (x->numInstances+1) * sizeof(t_instance));
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), (x->numInstances+1) * sizeof(int));
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, x->numInstances * sizeof(t_knn_info), (x->numInstances+1) * sizeof(t_knn_info));
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, x->numInstances * sizeof(int), (x->numInstances+1) * sizeof(int));
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numInstances * sizeof(t_member), (x->numInstances+1) * sizeof(t_member));
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numInstances * sizeof(int), (x->numInstances+1) * sizeof(int));
	
	x->instanceClusterMembership[instanceIdx] = instanceIdx;
	x->instanceFeatureLengths[instanceIdx] = listLength;
	x->clusterMemberLengths[instanceIdx] = 2; // 2 because we're unclustered to start, and each instance has a cluster with itself as a member, plus -1 as the 2nd element to terminate the list
    
	x->instances[instanceIdx].instance = (t_float *)t_getbytes_(listLength * sizeof(t_float));
    
	x->clusterMembers[instanceIdx].member = (int *)t_getbytes_(2 * sizeof(int));
	
	// init new clusterMembers
	x->clusterMembers[instanceIdx].member[0] = instanceIdx; // first member of the cluster is the instance index
	x->clusterMembers[instanceIdx].member[1] = -1;
	
	x->numInstances++;
	x->numClusters++;
	x->neighborhood++;
	
	if(x->featureLength != listLength)
	{
		x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));
		x->attributeOrder = (int *)t_resizebytes_(x->attributeOrder, x->featureLength * sizeof(int), listLength * sizeof(int));
		x->weights = (t_float *)t_resizebytes_(x->weights, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));
		x->featureLength = listLength;
        
		// initialize attributeOrder
		for(i=0; i<x->featureLength; i++)
			x->attributeOrder[i] = i;
        
		// initialize weights
		for(i=0; i<x->featureLength; i++)
			x->weights[i] = 1.0;
		
		x->attributelo = 0;
		x->attributehi = x->featureLength-1;
		post("feature length: %i.", x->featureLength);
		post("attribute range: %i through %i.", x->attributelo, x->attributehi);
	};
    
	for(i=0; i<listLength; i++)
		x->instances[instanceIdx].instance[i] = atom_getfloat(argv+i);
	
	outlet_int(x->id, instanceIdx); // output received feedback here rather than a post (in case of hi-speed training)
}


static void timbreID_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	t_float sum, best, second_best;
	t_float distanceOutput, confidence;
    int i, listLength, id, *votes, top_vote;
	//s=s; // to get rid of 'unused variable' warning
    
	if(x->numInstances)
	{
		votes = (int *)t_getbytes_(x->numClusters * sizeof(int));
		
		// init votes to 0
		for(i=0; i<x->numClusters; i++)
			votes[i] = 0;
        
		// init cluster info to instance idx
		for(i=0; i<x->numInstances; i++)
			x->knnDistsIdxs[i].cluster = i;
        
		listLength = argc;
		distanceOutput = 0;
		confidence = 0;
		
		if(x->featureLength != listLength)
		{
			x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));
			x->featureLength = listLength;
			x->attributelo = 0;
			x->attributehi = x->featureLength-1;
			post("feature length: %i.", x->featureLength);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
        
		for(i=0; i<x->featureLength; i++)
			x->featureInput[i] = atom_getfloat(argv+i);		
        
		id = 0;
		best = FLT_MAX;
		
		for(i=0; i<x->numInstances; i++)
		{
			sum = 0;		
            
			switch(x->distMetric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->featureInput, x->instances[i].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->featureInput, x->instances[i].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->featureInput, x->instances[i].instance);
					break;
				default:
					break;
			};
            
			x->knnDistsIdxs[i].dist = x->knnDistsIdxs[i].safeDist = sum; // store the distance
			x->knnDistsIdxs[i].idx = i; // store the idx
		};
        
		// a reduced sort, so that the first k elements in knnDistsIdxs will be the lowest distances in the list, and in order to boot.	
		timbreID_sort_knn_info(x->k, x->numInstances, -1, x->knnDistsIdxs); // pass a prevMatch value of -1, since it's unused here
        
		// store instance's cluster id
		for(i=0; i<x->k; i++)
			x->knnDistsIdxs[i].cluster = x->instanceClusterMembership[x->knnDistsIdxs[i].idx];
		
		// vote
		for(i=0; i<x->k; i++)
			votes[x->knnDistsIdxs[i].cluster]++;
        
		top_vote = -1;
		for(i=0; i<x->numClusters; i++)
			if(votes[i] > top_vote)
			{
				top_vote = votes[i];
				id = i; // store cluster id of winner
			};
		
		// in case of a tie, pick the shortest distance
		if(top_vote <= (x->k*0.5))
			id = x->knnDistsIdxs[0].cluster;
        
		for(i=0; i<x->k; i++)
			if(x->knnDistsIdxs[i].cluster==id)
			{
				best = x->knnDistsIdxs[i].safeDist;
				break;
			};
		
		second_best = FLT_MAX;
		
		for(i=0; i<x->k; i++)
			if(x->knnDistsIdxs[i].cluster!=id)
			{
				second_best = x->knnDistsIdxs[i].safeDist;
				break;
			};
        
		// if no second best assignment is made (because all K items belong to same cluster), make 2nd best the 2nd in list
		if(second_best==FLT_MAX)
			second_best = x->knnDistsIdxs[1].safeDist;
		
		if(second_best<=0 || second_best==FLT_MAX)
			confidence = 0;
		else
			confidence = 1-(best/second_best);
		
		distanceOutput = best;
        
		// free memory
		t_freebytes_(votes, x->numClusters*sizeof(int));
		
		outlet_float(x->confidence, confidence);	
		outlet_float(x->nearestDist, distanceOutput);
		outlet_int(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform ID.");
	
}


static void timbreID_worst_match(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	t_float sum, worst;
	t_float distanceOutput;
    int i, listLength, id;
	//s=s; // to get rid of 'unused variable' warning
    
	if(x->numInstances)
	{
		listLength = argc;
		distanceOutput = 0;
        
		if(x->featureLength != listLength)
		{
			x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));
			x->featureLength = listLength;
			x->attributelo = 0;
			x->attributehi = x->featureLength-1;
			post("feature length: %i.", x->featureLength);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
        
		for(i=0; i<argc; i++)
			x->featureInput[i] = atom_getfloat(argv+i);		
        
		id = 0;
		worst = 0;
		
		for(i=0; i<x->numInstances; i++)
		{
			sum = 0;		
            
			switch(x->distMetric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->featureInput, x->instances[i].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->featureInput, x->instances[i].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->featureInput, x->instances[i].instance);
					break;
				default:
					break;
			};
            
			if(sum > worst)
			{
				worst = sum;
				id = i;
				post("updated worst: %i, %f", id, worst);
			}
		};
        
		id = x->instanceClusterMembership[id];
		
		distanceOutput = worst;
		
		outlet_float(x->confidence, 0);	
		outlet_float(x->nearestDist, distanceOutput);
		outlet_int(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform worst match.");
	
}


static void timbreID_forget(t_timbreID *x)
{	
	if(x->numInstances > 0)
	{
		// free the instance
		t_freebytes_(x->instances[x->numInstances-1].instance, x->instanceFeatureLengths[x->numInstances-1]*sizeof(t_float));		
		x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), (x->numInstances-1) * sizeof(t_instance));
		
		// shrink instanceFeatureLengths
		x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), (x->numInstances-1) * sizeof(int));
		
		x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, x->numInstances * sizeof(t_knn_info), (x->numInstances-1) * sizeof(t_knn_info));
		x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, x->numInstances * sizeof(int), (x->numInstances-1) * sizeof(int));
        
		// should probably do a double-check here that the database isn't clustered. if it is clustered, we need to find
		// which cluster the instance belongs to and remove that instance from the cluster's member list.
		
		// shrink clusterMembers
		t_freebytes_(x->clusterMembers[x->numInstances-1].member, x->clusterMemberLengths[x->numInstances-1]*sizeof(int));		
		x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numInstances * sizeof(t_member), (x->numInstances-1) * sizeof(t_member));
		
		// shrink instanceFeatureLengths
		x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numInstances * sizeof(int), (x->numInstances-1) * sizeof(int));
        
		x->numInstances--;
		x->numClusters--;
        
		post("forgot last instance. instances 0 through %i remain.", x->numInstances-1);
	}
	else
		error("timbreID: nothing to forget.");
}



//************************************* concatenative synthesis functions

static void timbreID_concat_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	t_float sum, best;
	t_float distanceOutput;
    int i, j, searchStart, searchFinish, listLength, halfNeighborhood, id;
	//s=s; // to get rid of 'unused variable' warning
    
	if(x->numInstances)
	{
		// init knn info
		for(i=0; i<x->numInstances; i++)
		{
			x->knnDistsIdxs[i].idx = i;
			x->knnDistsIdxs[i].dist = x->knnDistsIdxs[i].safeDist = FLT_MAX;
		}
		
		listLength = argc;
		distanceOutput = 0;
		
		halfNeighborhood = x->neighborhood*0.5;
		
		if(x->featureLength != listLength)
		{
			x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), listLength * sizeof(t_float));
			x->featureLength = listLength;
			x->attributelo = 0;
			x->attributehi = x->featureLength-1;
			post("feature length: %i.", x->featureLength);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
        
		for(i=0; i<argc; i++)
			x->featureInput[i] = atom_getfloat(argv+i);		
        
		id = 0;
		best = FLT_MAX;
		
		// for just the searchStart, check to see if x->neighborhood was EVEN.  if so, we should make searchStart = searchCenter - halfNeighborhood + 1
		if(x->neighborhood%2==0)
			searchStart = x->searchCenter - halfNeighborhood + 1;
		else
			searchStart = x->searchCenter - halfNeighborhood;
		
		if(searchStart<0)
			searchStart = x->numInstances + searchStart;  // wraps in reverse to end of table.  + i because i is neg.
		else
			searchStart = searchStart%x->numInstances;
        
		searchFinish = x->searchCenter + halfNeighborhood;
		if(searchFinish<0)
			searchFinish = x->numInstances + searchFinish;  // wraps in reverse to end of table.  + i because i is neg.
		else
			searchFinish = searchFinish%x->numInstances;
        
		for(j=0, i=searchStart; j<x->neighborhood; j++)
		{
			sum = 0;
            
			switch(x->distMetric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->featureInput, x->instances[i].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->featureInput, x->instances[i].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->featureInput, x->instances[i].instance);
					break;
				default:
					break;
			};
            
            
			x->knnDistsIdxs[i].dist = x->knnDistsIdxs[i].safeDist = sum; // store the distance
            
			i++;
			i = i%x->numInstances;
		};
        
		// a reduced sort, so that the first maxMatches elements in knnDistsIdxs will be the lowest distances in the list, and in order to boot.
		// pass x->prevMatch to make sure we don't output the same match two times in a row (to prevent one grain being played back several
		// times in sequence.
		// this is wasteful in restricted searches because we don't need to look through all x->numInstances
		timbreID_sort_knn_info(x->maxMatches, x->numInstances, x->prevMatch, &x->knnDistsIdxs[0]);
        
		if(x->prevMatch == -1)
			id = x->knnDistsIdxs[0].idx;
		else
		{
			for(i=0, best=FLT_MAX; i<x->maxMatches; i++)
			{	
				t_float dist;
				dist=sum=0.0;
                
				switch(x->distMetric)
				{
					case 0:
						sum = timbreID_squared_euclid(x, x->instances[x->prevMatch].instance, x->instances[x->knnDistsIdxs[i].idx].instance);
						break;
					case 1:
						sum = timbreID_manhattan(x, x->instances[x->prevMatch].instance, x->instances[x->knnDistsIdxs[i].idx].instance);
						break;
					case 2:
						sum = timbreID_correlation(x, x->instances[x->prevMatch].instance, x->instances[x->knnDistsIdxs[i].idx].instance);
						break;
					default:
						break;
				};
                
				if( sum < best)
				{
					best = sum;
					id = x->knnDistsIdxs[i].idx;
				};
                
				sum = 0;
			};
		};
        
		if(x->reorientFlag)
			x->searchCenter = id;
		
		if( rand() < RAND_MAX*x->jumpProb )
		{
			id += rand();
			id = id%x->numInstances;
			x->searchCenter = id;
		};
        
		if(x->stutterProtect)
			x->prevMatch = id;
		else
			x->prevMatch = -1;
        
		distanceOutput = best;
        
		outlet_float(x->nearestDist, distanceOutput);
		outlet_int(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform ID.");
	
}


static void timbreID_concat_neighborhood(t_timbreID *x, int n)
{
	n = (n>x->numInstances)?x->numInstances:n;
	n = (n<1)?1:n;
	x->neighborhood = n;	
}


static void timbreID_concat_jumpProb(t_timbreID *x, int jp)
{
	jp = (jp<0)?0:jp;
	jp = (jp>1)?1:jp;
	x->jumpProb = jp;
}


static void timbreID_concat_reorient(t_timbreID *x, int r)
{
	r = (r<0)?0:r;
	r = (r>1)?1:r;
	x->reorientFlag = r;
}


static void timbreID_concat_searchCenter(t_timbreID *x, int sc)
{
	if( sc < 0 )
		x->searchCenter = 0;
	else if( sc >=x->numInstances )
		x->searchCenter = x->numInstances-1;
	else
		x->searchCenter = sc;
}


static void timbreID_concat_maxMatches(t_timbreID *x, int mm)
{	
	if (mm < 1) x->maxMatches = 1;
	else if (mm > 50) x->maxMatches = 50;
	else x->maxMatches = mm;
}


static void timbreID_concat_stutterProtect(t_timbreID *x, int sp)
{
	sp = (sp<0)?0:sp;
	sp = (sp>1)?1:sp;
	x->stutterProtect = sp;
	
	post("stutter protect: %i", x->stutterProtect);
}

//************************************* END concatenative synthesis functions




static void timbreID_knn(t_timbreID *x, int k)
{		
	if(k<1.0)
		post("k must be greater than zero.");
	else if(k>x->numInstances)
		post("k must be less than the total number of instances.");
	else
	{
		x->k = k;
		post("searching %i neighbors for KNN.", x->k);
	}
}


static void timbreID_normalize(t_timbreID *x, int n)
{	
	int i, j;
	t_float *attribute_column;
	
	// create local memory
	attribute_column = (t_float *)t_getbytes_(x->numInstances * sizeof(t_float));
	
	if(n<=0)
	{
		// free memory	
		t_freebytes_(x->normData, x->featureLength*sizeof(t_normData));
		
		x->normalize=0;
		post("feature attribute normalization OFF.");
	}
	else
	{
		if(x->numInstances)
		{
			// create memory;
			x->normData = (t_normData *)t_getbytes_(x->featureLength * sizeof(t_normData));
            
			// j for columns (attributes), i for rows (instances)
			for(j=0; j<x->featureLength; j++)
			{
				for(i=0; i<x->numInstances; i++)
					attribute_column[i] = x->instances[i].instance[j];
                
				timbreID_sort_float(x->numInstances, &attribute_column[0]);
				
				x->normData[j].min = attribute_column[0];
				x->normData[j].max = attribute_column[x->numInstances-1];
				
				// don't divide by zero
				if(x->normData[j].max <= x->normData[j].min)
				{
					x->normData[j].max = 2.0;
					x->normData[j].min = 1.0;
				};
				
				x->normData[j].denominator = 1.0/(x->normData[j].max - x->normData[j].min);
			};
			
			x->normalize=1;
			post("feature attribute normalization ON.");
		}
		else
			error("timbreID: no training instances have been loaded. cannot calculate normalization terms.");
	}
    
	// free local memory	
	t_freebytes_(attribute_column, x->numInstances*sizeof(t_float));
    
}


static void timbreID_print(t_timbreID *x)
{    
	post("no. of instances: %i", x->numInstances);
	post("feature length: %i", x->featureLength);
	post("attribute range: %i through %i\n", x->attributelo, x->attributehi);
    
	post("normalization: %i", x->normalize);
	post("distance metric: %i", x->distMetric);
	post("no. of clusters: %i", x->numClusters);
	post("KNN: %i\n", x->k);
	
	post("relative ordering: %i\n", x->relativeOrdering);
	
	post("search center: %i", x->searchCenter);
	post("neighborhood: %i", x->neighborhood);
	post("reorient: %i", x->reorientFlag);
	post("max matches: %i", x->maxMatches);
	post("jump probability: %i", x->jumpProb);
	post("stutter protect: %i", x->stutterProtect);
}


static void timbreID_numInstances(t_timbreID *x)
{
	t_atom *listOut;
    
	// create local memory
	listOut = (t_atom *)t_getbytes_(sizeof(t_atom));
    
	atom_setfloat(listOut, x->numInstances);
    
	outlet_list(x->x_featureList, 0, 1, listOut);
    
	// free local memory
	t_freebytes_(listOut, sizeof(t_atom));
}


static void timbreID_manual_cluster(t_timbreID *x, int numClusters, int cluster_idx, int low, int hi)
{
	int i, j, cluster_idx_i, low_i, hi_i, num_members;
    
	cluster_idx_i = cluster_idx;
	low_i = low;
	hi_i = hi;
	num_members = hi_i - low_i + 1;
    
	if(x->numInstances < numClusters)
		error("timbreID: not enough instances to cluster.");
	else
	{
		// only change memory size if x->numClusters hasn't been updated to be equal to numClusters
		if(x->numClusters != numClusters)
		{
			x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), numClusters * sizeof(t_member));
			x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), numClusters * sizeof(int));
            
			x->numClusters = numClusters;
		};
        
		// free the old memory for this cluster member list
		t_freebytes_(x->clusterMembers[cluster_idx_i].member, x->clusterMemberLengths[cluster_idx_i]*sizeof(int));
        
		// update the size of the list
		x->clusterMemberLengths[cluster_idx_i] = num_members+1; // +1 for the terminating -1
        
		x->clusterMembers[cluster_idx_i].member = (int *)t_getbytes_(x->clusterMemberLengths[cluster_idx_i] * sizeof(int));		
        
		for(i=low_i; i<=hi_i; i++)
			x->instanceClusterMembership[i] = cluster_idx_i;
        
		for(i=low_i, j=0; i<=hi_i; i++, j++)
			x->clusterMembers[cluster_idx_i].member[j] = i;
        
		// terminate with -1
		x->clusterMembers[cluster_idx_i].member[j] = -1;
        
		post("cluster %i contains instances %i through %i.", cluster_idx_i, low_i, hi_i);
	};
}


static void timbreID_compute_cluster(t_timbreID *x, int numClusters)
{
	int i, j, k, numInstances, numInstancesM1, num_pairs, num_clusterMembers1, num_clusterMembers2, num_clusterMembers_sum, clusterCount, *min_dist_idx;
	t_instance *cluster_data;
	t_float *pair_dists;
	t_float min_dist, num_clusterMembers1_recip;
	t_atom *listOut;
    
	if(x->numInstances < numClusters)
		error("timbreID: not enough instances to cluster.");
	else if(x->numClusters != x->numInstances)
		error("timbreID: instances already clustered. uncluster first.");
	else if(numClusters == 0)
		error("timbreID: cannot create 0 clusters.");
	else
	{
        x->numClusters = numClusters;
        numInstances = x->numInstances;
        numInstancesM1 = numInstances-1;
        num_pairs = (numInstances*numInstancesM1) * 0.5;
        clusterCount = numInstances;
        num_clusterMembers1 = 0;
        num_clusterMembers2 = 0;
        num_clusterMembers1_recip = 1;
        i=j=k=0;
        
        // create local memory
        min_dist_idx = (int *)t_getbytes_(2 * sizeof(int));
        cluster_data = (t_instance *)t_getbytes_(numInstances * sizeof(t_instance));
        pair_dists = (t_float *)t_getbytes_(num_pairs * sizeof(t_float));
        listOut = (t_atom *)t_getbytes_(numInstances * sizeof(t_atom));
        
        for(i=0; i<numInstances; i++)
        {	
            x->clusterMembers[i].member[0] = i; // first member of the cluster is the instance index
            x->clusterMembers[i].member[1] = -1;
        }
        
        // copy x->instances into a safe local copy: cluster_data
        for(i=0; i<numInstances; i++)
        {
            cluster_data[i].instance = (t_float *)t_getbytes_(x->featureLength * sizeof(t_float));
            
            for(j=0; j<x->featureLength; j++)
                cluster_data[i].instance[j] = x->instances[i].instance[j];
        }
        
        
        while(clusterCount > x->numClusters)
        {		
            min_dist = FLT_MAX;
            
            // init min_dist_idx
            for(i=0; i<2; i++)
                min_dist_idx[i] = -1;
            
            // init pair distances 
            for(i=0; i<num_pairs; i++)
                pair_dists[i] = FLT_MAX;
			
            
            // get distances between all possible pairs in cluster_data
            for(i=0, k=0; i<numInstancesM1; i++)
            {
                if( cluster_data[i].instance[0] != -9999 ) // if this is true, the data hasn't been clustered yet.
                {
                    for(j=1; j<numInstances; j++)
                    {	
                        if( (i+j) < numInstances )
                        {
                            if( cluster_data[i+j].instance[0] != -9999 )
                            {
                                switch(x->distMetric)
                                {
                                    case 0:
                                        pair_dists[k] = timbreID_squared_euclid(x, cluster_data[i].instance, cluster_data[i+j].instance);
                                        break;
                                    case 1:
                                        pair_dists[k] = timbreID_manhattan(x, cluster_data[i].instance, cluster_data[i+j].instance);
                                        break;
                                    case 2:
                                        pair_dists[k] = timbreID_correlation(x, cluster_data[i].instance, cluster_data[i+j].instance);
                                        break;
                                    default:
                                        break;
                                };
                                
                                num_clusterMembers1 = x->clusterMemberLengths[i]-1; // -1 because the list is terminated with -1
                                num_clusterMembers2 = x->clusterMemberLengths[i+j]-1;
								
                                // definition of Ward's linkage from MATLAB linkage doc
                                // pair_dists[k] is already squared euclidean distance
                                
                                num_clusterMembers_sum = num_clusterMembers1 + num_clusterMembers2;
                                
                                if(num_clusterMembers_sum > 0)
                                    pair_dists[k] = num_clusterMembers1*num_clusterMembers2 * (pair_dists[k]/(num_clusterMembers1+num_clusterMembers2));
                                else
                                    pair_dists[k] = FLT_MAX;
                                
                                if(pair_dists[k]<min_dist)
                                {
                                    min_dist=pair_dists[k];
                                    min_dist_idx[0]=i;
                                    min_dist_idx[1]=i+j;
                                };
                                
                                k++; // increment pair_dists index if something was actually written to it.
                            };
                        }
                        else
                            break;
                    }
                }
            };
            
            // we've found the smallest distance between cluster_data elements and stored it 
            // in min_dist. we've store the cluster_data indices of the two elements in 
            // min_dist_idx[0] and min_dist_idx[1].
            
            // set i to the index for storing the new member(s) of the cluster.
            i = x->clusterMemberLengths[min_dist_idx[0]]-1;
            
            // actually store the new member(s).
            j=0;
            while(x->clusterMembers[min_dist_idx[1]].member[j] != -1)
            {
                // make some more memory for the new member(s)
                x->clusterMembers[min_dist_idx[0]].member = (int *)t_resizebytes_(x->clusterMembers[min_dist_idx[0]].member, x->clusterMemberLengths[min_dist_idx[0]] * sizeof(int), (x->clusterMemberLengths[min_dist_idx[0]]+1) * sizeof(int));
                x->clusterMemberLengths[min_dist_idx[0]]++; // remember to update this member list's length
                
                x->clusterMembers[min_dist_idx[0]].member[i++] = x->clusterMembers[min_dist_idx[1]].member[j++];
            }
            
            i = x->clusterMemberLengths[min_dist_idx[0]]-1;
            x->clusterMembers[min_dist_idx[0]].member[i] = -1; // terminate
            
            num_clusterMembers1 = x->clusterMemberLengths[min_dist_idx[0]]-1;
            
            if(num_clusterMembers1 > 0)
                num_clusterMembers1_recip = 1.0/(t_float)num_clusterMembers1;
            else
                num_clusterMembers1_recip = 1.0;
            
            // resize the usurped cluster's cluster list memory, and update its size to 1
            x->clusterMembers[min_dist_idx[1]].member = (int *)t_resizebytes_(x->clusterMembers[min_dist_idx[1]].member, x->clusterMemberLengths[min_dist_idx[1]] * sizeof(int), sizeof(int));
            x->clusterMembers[min_dist_idx[1]].member[0] = -1;
            x->clusterMemberLengths[min_dist_idx[1]] = 1;
            
            // grab the first original instance for this cluster index
            for(i=0; i<x->featureLength; i++)
                cluster_data[min_dist_idx[0]].instance[i] = x->instances[min_dist_idx[0]].instance[i];
            
            // sum the original instances of the cluster members to compute centroid below
            for(i=1; i<num_clusterMembers1; i++)
                for(j=0; j<x->featureLength; j++)
                    cluster_data[min_dist_idx[0]].instance[j] += x->instances[  x->clusterMembers[min_dist_idx[0]].member[i]  ].instance[j];
            
            // compute centroid
            for(i=0; i<x->featureLength; i++)
                cluster_data[min_dist_idx[0]].instance[i] *= num_clusterMembers1_recip;
            
            // write -9999 to the first element in the nearest neighbor's instance to indicate it's now vacant.
            // this is all that's needed since all previous members were averaged and stored here.
            cluster_data[min_dist_idx[1]].instance[0] = -9999;
            
            clusterCount--;
        };
        
        // since the indices of the clusters have gaps from the process,
        // shift the clusterMembers arrays that actually have content (!= -1)
        // to the head of clusterMembers.  this will produce indices from 0 through numClusters-1.
        for(i=0, k=0; i<numInstances; i++)
            if( x->clusterMembers[i].member[0] != -1)
            {
                // resize this member list
                x->clusterMembers[k].member = (int *)t_resizebytes_(x->clusterMembers[k].member, x->clusterMemberLengths[k] * sizeof(int), x->clusterMemberLengths[i] * sizeof(int));
                
                for(j=0; j<x->clusterMemberLengths[i]; j++)
                    x->clusterMembers[k].member[j] = x->clusterMembers[i].member[j];
                
                // shift the list length info back
                x->clusterMemberLengths[k] = x->clusterMemberLengths[i];
                
                k++;
            };
        
        // free the excess clusterMembers memory
        for(i=x->numClusters; i<numInstances; i++)
            t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
		
        // resize clusterMembers so it is only x->numClusters big
        x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, numInstances * sizeof(t_member), x->numClusters * sizeof(t_member));
        
        // resize clusterMemberLengths so it is only x->numClusters big
        x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, numInstances * sizeof(int), x->numClusters * sizeof(int));
		
        for(i=0, k=0; i<x->numClusters; i++)
            for(j=0; j<(x->clusterMemberLengths[i]-1); j++)
            {
                x->instanceClusterMembership[x->clusterMembers[i].member[j]] = i;
                atom_setfloat(listOut+k, x->clusterMembers[i].member[j]);
                k++;
            };
		
        outlet_list(x->x_orderList, 0, x->numInstances, listOut);
        
        // free memory
        t_freebytes_(min_dist_idx, 2*sizeof(int));
        
        // free the database memory
        for(i=0; i<numInstances; i++)
            t_freebytes_(cluster_data[i].instance, x->featureLength*sizeof(t_float));
		
        t_freebytes_(cluster_data, numInstances*sizeof(t_instance));
        
        t_freebytes_(pair_dists, num_pairs*sizeof(t_float));
        t_freebytes_(listOut, numInstances*sizeof(t_atom));
		
        post("instances clustered.");
        
	} // end of main if/else	
}


static void timbreID_uncluster(t_timbreID *x)
{
	int i;
    
	// free each x->clusterMembers list's memory
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), x->numInstances * sizeof(t_member));
    
	for(i=0; i<x->numInstances; i++)
		x->clusterMembers[i].member = (int *)t_getbytes_(2 * sizeof(int));
    
	// expand size of clusterMemberLengths again
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), x->numInstances * sizeof(int));
    
	x->numClusters = x->numInstances;
	
	for(i=0; i<x->numInstances; i++)
	{
		x->instanceClusterMembership[i]=i; // init membership to index
		x->clusterMembers[i].member[0]=i; // first member of the cluster is the instance index
		x->clusterMembers[i].member[1] = -1;
        
		x->clusterMemberLengths[i] = 2;
	}
	
    post("instances unclustered.");
}


static void timbreID_compute_variance(t_timbreID *x)
{
	int i, j;
	t_float max, *attribute_var;
	t_instance *meanCentered;
	
	if(x->numInstances > 0)
	{
		// create local memory
		attribute_var = (t_float *)t_getbytes_(x->featureLength * sizeof(t_float));
		meanCentered = (t_instance *)t_getbytes_(x->numInstances * sizeof(t_instance));
        
		for(i=0; i<x->numInstances; i++)
			meanCentered[i].instance = (t_float *)t_getbytes_(x->featureLength * sizeof(t_float));
		
		// init mean centered
		for(i=0; i<x->numInstances; i++)
			for(j=0; j<x->featureLength; j++)
				meanCentered[i].instance[j] = 0.0;
		
		// get the mean of each attribute
		// mean() checks for FLT_MAX and doesn't include those rows
		for(i=0; i<x->featureLength; i++)
			attribute_var[i] = timbreID_mean(x->numInstances, i, x->instances, x->normalize, x->normData);
		
		// center the data and write the matrix B
		for(i=0; i<x->numInstances; i++)
			for(j=0; j<x->featureLength; j++)
			{
				if(x->normalize)
					meanCentered[i].instance[j] = ((x->instances[i].instance[j] - x->normData[j].min) * x->normData[j].denominator) - attribute_var[j];
				else
					meanCentered[i].instance[j] = x->instances[i].instance[j] - attribute_var[j];				
			}
        
        // 	// variance is calculated as: sum(B(:,1).^2)/(M-1) for the first attribute
        // 	// run process by matrix columns rather than rows, hence the j, i order
		for(j=0; j<x->featureLength; j++)
		{		
			attribute_var[j] = 0;
			
			for(i=0; i<x->numInstances; i++)
				if(x->instances[i].instance[0] == FLT_MAX)
					continue;
				else
					attribute_var[j] += meanCentered[i].instance[j] * meanCentered[i].instance[j];
			
			if((x->numInstances-1) > 0)
				attribute_var[j] /= x->numInstances-1;
		}
        
        // 	for(i=0; i<x->featureLength; i++)
        // 		post("attribute variance %i: %f.", i, attribute_var[i]);
		
		// sort attributeOrder by largest variances: find max in attribute_var,
		// replace it with -99999, find next max.
		for(i=0; i<x->featureLength; i++)
		{
			max=0.0;
			for(j=0; j<x->featureLength; j++)
			{   
				if(attribute_var[j] > max)
				{
					max = attribute_var[j];
					x->attributeOrder[i] = j;
				}
			};
			
			attribute_var[x->attributeOrder[i]] = -99999.0;
		};
		
		// free local memory
		t_freebytes_(attribute_var, x->featureLength*sizeof(t_float));
        
		// free the meanCentered memory
		for(i=0; i<x->numInstances; i++)
			t_freebytes_(meanCentered[i].instance, x->featureLength*sizeof(t_float));
		
		t_freebytes_(meanCentered, x->numInstances*sizeof(t_instance));
        
		post("attributes ordered by variance.");
	}
	else
		error("timreID: no instances for variance computation.");
    
}


static void timbreID_clusters_list(t_timbreID *x)
{
	int i, j, k;
	t_atom *listOut;
    
	// create local memory
	listOut = (t_atom *)t_getbytes_(x->numInstances * sizeof(t_atom));
    
	for(i=0, k=0; i<x->numClusters; i++)
		for(j=0; j<(x->clusterMemberLengths[i]-1); j++, k++) // -1 because it's terminated by -1
			atom_setfloat(listOut+k, x->clusterMembers[i].member[j]);
    
	outlet_list(x->x_orderList, 0, x->numInstances, listOut);
	
	// free local memory
	t_freebytes_(listOut, x->numInstances*sizeof(t_atom));
}


static void timbreID_cluster_list(t_timbreID *x, int idx)
{
	int i, idx_i;
	t_atom *listOut;
    
	idx_i = idx;
    
	if(idx_i >= x->numClusters || idx_i < 0)
		error("timbreID: cluster %i does not exist.", idx_i);
	else
	{
		// create local memory
		listOut = (t_atom *)t_getbytes_(0);
        
		for(i=0; i<(x->clusterMemberLengths[idx_i]-1); i++)
		{
			listOut = (t_atom *)t_resizebytes_(listOut, i * sizeof(t_atom), (i+1) * sizeof(t_atom));
			atom_setfloat(listOut+i, x->clusterMembers[idx_i].member[i]);
		};
        
		outlet_list(x->x_orderList, 0, i, listOut);
		
		// free local memory
		t_freebytes_(listOut, i*sizeof(t_atom));
	}
}


static void timbreID_clusterMembership(t_timbreID *x, int idx)
{
	int idx_i;
	t_atom *listOut;
    
	idx_i = idx;
    
	if(idx_i >= x->numInstances || idx_i < 0)
		error("timbreID: instance %i does not exist.", idx_i);
	else
	{
		// create local memory for a single element
		listOut = (t_atom *)t_getbytes_(sizeof(t_atom));
        
		atom_setfloat(listOut, x->instanceClusterMembership[idx_i]);
		outlet_list(x->x_orderList, 0, 1, listOut);
		
		// free local memory
		t_freebytes_(listOut, sizeof(t_atom));
	}
}


static void timbreID_compute_order(t_timbreID *x, int reference)
{
	int i, j, smallIdx, ref;
	t_float smallest, sum;
	t_instance *instances;
	t_atom *listOut;
    
	// create local memory
	instances = (t_instance *)t_getbytes_(x->numInstances * sizeof(t_instance));
	listOut = (t_atom *)t_getbytes_(x->numInstances * sizeof(t_atom));
    
	for(i=0; i<x->numInstances; i++)
		instances[i].instance = (t_float *)t_getbytes_(x->featureLength * sizeof(t_float));
	
	if(reference >= x->numInstances)
		ref = x->numInstances-1;
	else if(reference < 0)
		ref = 0;
	else
	    ref = reference;
    
	// make a local copy of instances so they can be abused
    for(i=0; i<x->numInstances; i++)
    	for(j=0; j<x->featureLength; j++)
    		instances[i].instance[j] = x->instances[i].instance[j];
    
    
    for(i=0; i<x->numInstances; i++)
    {
		smallest = FLT_MAX;
		smallIdx = 0;		
		
		for(j=0; j<x->numInstances; j++)
		{
			sum = 0;	
            
			// break out of this for iteration early if this instance slot has already been used.
			if(instances[j].instance[0] == FLT_MAX)
				continue;
            
			switch(x->distMetric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->instances[ref].instance, instances[j].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->instances[ref].instance, instances[j].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->instances[ref].instance, instances[j].instance);
					break;
				default:
					break;
			};
            
			if(sum<smallest)
			{
				smallest = sum;
				smallIdx = j;
			};
			
		};
        
		atom_setfloat(listOut+i, smallIdx); // store the best from this round;
        
		if(x->relativeOrdering)
			ref = smallIdx; // reorient search to nearest match;
		
		// set this instance to something huge so it will never be chosen as a good match
		for(j=0; j<x->featureLength; j++)
			instances[smallIdx].instance[j] = FLT_MAX;
        
	};
    
	outlet_list(x->x_orderList, 0, x->numInstances, listOut);
    
	// free local memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(instances[i].instance, x->featureLength*sizeof(t_float));
    
	t_freebytes_(instances, x->numInstances*sizeof(t_instance));
	t_freebytes_(listOut, x->numInstances*sizeof(t_atom));
}


static void timbreID_relativeOrdering(t_timbreID *x, int rel)
{
	if(rel<0)
		x->relativeOrdering = 0;
	else if (rel>1)
		x->relativeOrdering = 1;
	else
		x->relativeOrdering = rel;
	
	if(x->relativeOrdering)
		post("relative ordering ON.");
	else
		post("relative ordering OFF.");
}


static void timbreID_distMetric(t_timbreID *x, int f)
{		
	if(f < 0)
		x->distMetric = 0;
	if(f > 3)
		x->distMetric = 3;
	else
		x->distMetric = f;
    
	switch(x->distMetric)
	{
		case 0:
			post("distance metric: EUCLIDEAN.");
			break;
		case 1:
			post("distance metric: MANHATTAN (taxicab distance).");
			break;
		case 2:
			post("distance metric: PEARSON CORRELATION COEFF.");
			break;
		default:
			break;
	};
}



static void timbreID_weights(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i;
	//s=s; // to get rid of 'unused variable' warning
	
	if(argc > x->featureLength)
	{
		post("WARNING: weights list longer than current feature length");
		argc = x->featureLength;
	}
	
	for(i=0; i<argc; i++)
		x->weights[i] = atom_getfloat(argv+i);
    
	// if only the first few of a long feature vector are specified, fill in the rest with 1.0
	for(i=argc; i<x->featureLength; i++)
		x->weights[i] = 1.0;
}


static void timbreID_attributes(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i;
	//s=s; // to get rid of 'unused variable' warning
	
	if(argc > x->featureLength)
	{
		post("WARNING: attribute list longer than timbreID's current feature length");
		argc = x->featureLength;
	}
	else
		post("attribute list received.");
    
	
	for(i=0; i<argc; i++)
		x->attributeOrder[i] = atom_getfloat(argv+i);
    
	// fill any remainder with attribute 0
	for(i=argc; i<x->featureLength; i++)
		x->attributeOrder[i] = 0;
}


static void timbreID_attribute_range(t_timbreID *x, int lo, int hi)
{
    
	if(lo < x->featureLength)
		x->attributelo = lo;
	else
		x->attributelo = x->featureLength-1;
	
	
	if(hi < x->featureLength)
		x->attributehi = hi;
	else
		x->attributehi = x->featureLength-1;
	
	
	if(lo > hi)
		post("WARNING: low attribute > high attribute.  Correct this before further use.");
    
    post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    
}


static void timbreID_order_attributes(t_timbreID *x)
{
	int i;
    
	// initialize attributeOrder
	for(i=0; i<x->featureLength; i++)
		x->attributeOrder[i] = i;
    
	post("attribute order initialized.");
}


static void timbreID_print_instance(t_timbreID *x, int e, int f, int g)
{
	int i;
    
    
	if( e >= x->numInstances || e < 0)
	{
		error("instance %i does not exist", e);
	}
	else
	{
		post("T%i = [", e);
        
		for(i=f; i<(g+1); i++)
		{
			
			if(i != g)
				post("%f, ", x->instances[e].instance[i]);
			else
				post("%f", x->instances[e].instance[i]);
		};
        
        
		post("]");
		post("");
    };
    
}


static void timbreID_feature_list(t_timbreID *x, int idx)
{
	int i, idx_i;
	t_atom *listOut;
    
	idx_i = idx;
    
	if(idx_i >= x->numInstances || idx_i < 0)
		error("timbreID: instance %i does not exist.", idx_i);
	else
	{
		// create local memory
		listOut = (t_atom *)t_getbytes_(x->featureLength * sizeof(t_atom));
		
		for(i=0; i<x->featureLength; i++)
		{
			if(x->normalize)
			{
				if( x->normData[i].max <= x->normData[i].min )
				{
					x->normData[i].max = 2.0;
					x->normData[i].min = 1.0;
				}
				
				atom_setfloat(listOut+i, (x->instances[idx_i].instance[i] - x->normData[i].min)/(x->normData[i].max - x->normData[i].min));	
			}
			else
				atom_setfloat(listOut+i, x->instances[idx_i].instance[i]);
		}
		
		outlet_list(x->x_featureList, 0, x->featureLength, listOut);
		
		// free local memory
		t_freebytes_(listOut, x->featureLength*sizeof(t_atom));
	}
}


static void timbreID_similarityMatrix(t_timbreID *x, int startInstance, int finishInstance, int normalize)
{
	int i, j, k, l, numInst, startInst, finishInst, norm;
    
	if(x->numInstances)
	{
        startInstance = (startInstance<0)?0:startInstance;
        startInstance = (startInstance>=x->numInstances)?x->numInstances-1:startInstance;
        startInst = startInstance;
        
        finishInstance = (finishInstance<0)?0:finishInstance;
        finishInstance = (finishInstance>=x->numInstances)?x->numInstances-1:finishInstance;
        finishInst = finishInstance;
        
        normalize = (normalize<0)?0:normalize;
        normalize = (normalize>1)?1:normalize;
        norm = normalize;
        
        if(startInst>finishInst)
        {
            int tmp;
            
            tmp = finishInst;
            finishInst = startInst;
            startInst = tmp;
        }
        
        numInst = finishInst-startInst+1;
        
        if(numInst)
        {
            t_float maxDist;
            t_instance *distances;
            
            // create local memory
            distances = (t_instance *)t_getbytes_(numInst*sizeof(t_instance));
            
            for(i=0; i<numInst; i++)
            {
                distances[i].instance = (t_float *)t_getbytes_(numInst*sizeof(t_float));
                for(j=0; j<numInst; j++)
                    distances[i].instance[j] = 0.0;
            }
            
            maxDist = -1;
            
            for(i=startInst, j=0; i<=finishInst; i++, j++)
            {
                for(k=startInst, l=0; k<=finishInst; k++, l++)
                {
                    t_float dist;
                    
                    dist = 0;
                    
                    switch(x->distMetric)
                    {
                        case 0:
                            dist = timbreID_squared_euclid(x, x->instances[i].instance, x->instances[k].instance);
                            break;
                        case 1:
                            dist = timbreID_manhattan(x, x->instances[i].instance, x->instances[k].instance);
                            break;
                        case 2:
                            dist = timbreID_correlation(x, x->instances[i].instance, x->instances[k].instance);
                            break;
                        default:
                            break;
                    };
                    
                    if(dist>maxDist)
                        maxDist = dist;
                    
                    distances[j].instance[l] = dist;
                }
            }
            
            maxDist = 1.0/maxDist;
            
            for(i=startInst; i<numInst; i++)
            {
                t_atom *listOut;
                
                listOut = (t_atom *)t_getbytes_(numInst*sizeof(t_atom));
                
                for(j=0; j<numInst; j++)
                {
                    t_float dist;
                    
                    dist = distances[i].instance[j];
                    
                    if(norm)
                        dist *= maxDist;
                    
                    atom_setfloat(listOut+j, dist);
                }
                
                outlet_list(x->x_featureList, 0, numInst, listOut);
                
                t_freebytes_(listOut, numInst*sizeof(t_atom));
            }
            
            
            // free local memory
            for(i=0; i<numInst; i++)
                t_freebytes_(distances[i].instance, numInst*sizeof(t_float));
			
            t_freebytes_(distances, numInst*sizeof(t_instance));
            
        }
        else
            error("timbreID: bad range of instances");
	}
	else
		error("timbreID: no training instances have been loaded.");
}


static void timbreID_max_values(t_timbreID *x)
{
	int i;
	t_atom *listOut;
    
	if(x->normalize)
	{
		// create local memory
		listOut = (t_atom *)t_getbytes_(x->featureLength * sizeof(t_atom));
		
		for(i=0; i<x->featureLength; i++)
			atom_setfloat(listOut+i, x->normData[i].max);
		
		outlet_list(x->x_featureList, 0, x->featureLength, listOut);
		
		// free local memory
		t_freebytes_(listOut, x->featureLength*sizeof(t_atom));
	}
	else
		error("timbreID: feature database not normalized yet");
}


static void timbreID_min_values(t_timbreID *x)
{
	int i;
	t_atom *listOut;
    
	if(x->normalize)
	{
		// create local memory
		listOut = (t_atom *)t_getbytes_(x->featureLength * sizeof(t_atom));
		
		for(i=0; i<x->featureLength; i++)
			atom_setfloat(listOut+i, x->normData[i].min);
		
		outlet_list(x->x_featureList, 0, x->featureLength, listOut);
		
		// free local memory
		t_freebytes_(listOut, x->featureLength*sizeof(t_atom));
	}
	else
		error("timbreID: feature database not normalized yet");
}


static void timbreID_clear(t_timbreID *x)
{
	int i;
	
	// free the database memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
	
	x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), 0);
    
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), 0);
	
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, x->numInstances * sizeof(t_knn_info), 0);
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, x->numInstances * sizeof(int), 0);
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), 0);
	x->numInstances = 0;
	x->neighborhood = 0;
	x->numClusters = 0;
	
    post("all instances cleared.");
}


static void timbreID_write(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *header;
    t_float *fp;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    
	header = (int *)t_getbytes_((x->numInstances+1) * sizeof(int)); // record the size of each instance's feature, plus 1 for numInstances
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "wb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	header[0] = x->numInstances;
	for(i=0; i<x->numInstances; i++)
		header[i+1] = x->instanceFeatureLengths[i];
    
	fwrite(header, sizeof(int), x->numInstances+1, fd);
    
    for(i=0; i<x->numInstances; i++)
    {    
		fp = x->instances[i].instance;
		fwrite(fp, sizeof(t_float), x->instanceFeatureLengths[i], fd);
   	};
   	
   	
    post("wrote %i instances to file: %s.", x->numInstances, buf);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
    t_freebytes_(header, (x->numInstances+1)*sizeof(int));
}


static void timbreID_read(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    t_float *fp;
	int i;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
    // erase old instances & clusters and resize to 0.
    
	// free the database memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
    
    x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), 0);
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), 0);
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, x->numInstances * sizeof(t_knn_info), 0);
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, x->numInstances * sizeof(int), 0);
    
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), 0);
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), 0);
    
	x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), 0);   
	x->attributeOrder = (int *)t_resizebytes_(x->attributeOrder, x->featureLength * sizeof(int), 0);
	x->weights = (t_float *)t_resizebytes_(x->weights, x->featureLength * sizeof(t_float), 0);
    
    fd = fopen(buf, "rb");
    
    if (!fd)
    {
        post("%s: open failed", buf);
        return;
    }
    
	fread(&x->numInstances, sizeof(int), 1, fd);
	
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, 0, x->numInstances * sizeof(int));
    
	fread(x->instanceFeatureLengths, sizeof(int), x->numInstances, fd);
	
	// should search for the min and max instance sizes and store them both.
	// could resize relevant memory according to max, but read point limits according to min.
	// for now, just assume they're all the same.
	x->featureLength = x->instanceFeatureLengths[0];	
	x->neighborhood = x->numInstances;
	x->numClusters = x->numInstances;
	x->attributelo = 0;
	x->attributehi = x->featureLength-1;
    
    // resize instances & clusterMembers to numInstances
    x->instances = (t_instance *)t_resizebytes_(x->instances, 0, x->numInstances*sizeof(t_instance));
    
	for(i=0; i<x->numInstances; i++)
		x->instances[i].instance = (t_float *)t_getbytes_(x->instanceFeatureLengths[i] * sizeof(t_float));
	
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, 0, x->numInstances * sizeof(t_knn_info));
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, 0, x->numInstances * sizeof(int));
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, 0, x->numInstances * sizeof(t_member));
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, 0, x->numInstances * sizeof(int));
	x->featureInput = (t_float *)t_resizebytes_(x->featureInput, 0, x->featureLength * sizeof(t_float));
	x->attributeOrder = (int *)t_resizebytes_(x->attributeOrder, 0, x->featureLength * sizeof(int));
	x->weights = (t_float *)t_resizebytes_(x->weights, 0, x->featureLength * sizeof(t_float));
    
	// initialize attributeOrder
	for(i=0; i<x->featureLength; i++)
		x->attributeOrder[i] = i;
    
	// initialize weights
	for(i=0; i<x->featureLength; i++)
		x->weights[i] = 1.0;
    
	// initialize feature input buffer
	for(i=0; i<x->featureLength; i++)
		x->featureInput[i] = 0.0;
	
	for(i=0; i<x->numInstances; i++)
	{
		x->clusterMembers[i].member = (int *)t_getbytes_(2 * sizeof(int));
        
		x->instanceClusterMembership[i]=i; // init membership to index
		x->clusterMembers[i].member[0]=i; // first member of the cluster is the instance index
		x->clusterMembers[i].member[1] = -1;
        
		x->clusterMemberLengths[i] = 2;
	}
    
    
	// finally, read in the instance data
    for(i=0; i<x->numInstances; i++)
    {
		fp = x->instances[i].instance;
		fread(fp, sizeof(t_float), x->instanceFeatureLengths[i], fd);
    };
    
    
    post("read %i instances from file: %s.\n", x->numInstances, buf);    
    post("feature length: %i.", x->featureLength);
	post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    int i, j, *header;
    t_float *fp;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
	header = (int *)t_getbytes_(2 * sizeof(int));
	
	j=0; // to keep track of no. of instances written.
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	// unlike the binary write/read, here we just assume a common feature length, and never look at instanceFeatureLengths...for now.
	header[0] = x->numInstances;
	header[1] = x->featureLength;
    
	for(i=0; i<2; i++)
		fprintf(fd, "%i ", header[i]);
    
	fprintf(fd, "\n\n");
    
    for(i=0; i<x->numInstances; i++)
    {    
		fp = x->instances[i].instance;
        
		j=0;
        
		// only write actual values, not the FLT_MAX placeholders
		while(j<x->instanceFeatureLengths[i])
		{
			fprintf(fd, "%6.20f ", *fp++);
			j++;
		};
		
		fprintf(fd, "\n\n");
   	};
    
    post("wrote %i instances to file: %s.", x->numInstances, buf);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
    t_freebytes_(header, 2*sizeof(int));
    
    fclose(fd);
}


static void timbreID_read_text(t_timbreID *x, t_symbol *s)
{
    
    FILE *fd;
    
    int i, j, *header;
    
    t_float *fp;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
	header = (int *)t_getbytes_(2 * sizeof(int));
	
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	// free the database memory
	for(i=0; i<x->numInstances; i++)
		t_freebytes_(x->instances[i].instance, x->instanceFeatureLengths[i]*sizeof(t_float));
    
    x->instances = (t_instance *)t_resizebytes_(x->instances, x->numInstances * sizeof(t_instance), 0);
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, x->numInstances * sizeof(int), 0);
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, x->numInstances * sizeof(t_knn_info), 0);
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, x->numInstances * sizeof(int), 0);
    
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), 0);
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), 0);
    
	x->featureInput = (t_float *)t_resizebytes_(x->featureInput, x->featureLength * sizeof(t_float), 0);   
	x->attributeOrder = (int *)t_resizebytes_(x->attributeOrder, x->featureLength * sizeof(int), 0);
	x->weights = (t_float *)t_resizebytes_(x->weights, x->featureLength * sizeof(t_float), 0);
    
    fd = fopen(buf, "r");
    
    if (!fd)
    {
        post("%s: open failed", buf);
        return;
    }
    
	//// unlike the binary write/read, here we just assume a common feature length, and never look at instanceFeatureLengths...for now.
	for(i=0; i<2; i++, header++)
		fscanf(fd, "%i", header);
	
	// in reverse order due to ptr arithmetic
	x->featureLength = *(--header);
	x->numInstances = *(--header);
    
	x->neighborhood = x->numInstances;
	x->numClusters = x->numInstances;
	x->attributelo = 0;
	x->attributehi = x->featureLength-1;
    
    // resize instances & clusterMembers to numInstances
    x->instances = (t_instance *)t_resizebytes_(x->instances, 0, x->numInstances*sizeof(t_instance));
    
	for(i=0; i<x->numInstances; i++)
		x->instances[i].instance = (t_float *)t_getbytes_(x->featureLength * sizeof(t_float));
    
	x->instanceFeatureLengths = (int *)t_resizebytes_(x->instanceFeatureLengths, 0, x->numInstances * sizeof(int));
    
	x->knnDistsIdxs = (t_knn_info *)t_resizebytes_(x->knnDistsIdxs, 0, x->numInstances * sizeof(t_knn_info));
	x->instanceClusterMembership = (int *)t_resizebytes_(x->instanceClusterMembership, 0, x->numInstances * sizeof(int));
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, 0, x->numInstances * sizeof(t_member));
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, 0, x->numInstances * sizeof(int));
	x->featureInput = (t_float *)t_resizebytes_(x->featureInput, 0, x->featureLength * sizeof(t_float));
	x->attributeOrder = (int *)t_resizebytes_(x->attributeOrder, 0, x->featureLength * sizeof(int));
	x->weights = (t_float *)t_resizebytes_(x->weights, 0, x->featureLength * sizeof(t_float));
    
	// initialize instance sizes
	for(i=0; i<x->numInstances; i++)
		x->instanceFeatureLengths[i] = x->featureLength;
    
	// initialize attributeOrder
	for(i=0; i<x->featureLength; i++)
		x->attributeOrder[i] = i;
    
	// initialize weights
	for(i=0; i<x->featureLength; i++)
		x->weights[i] = 1.0;
    
	// initialize feature input buffer
	for(i=0; i<x->featureLength; i++)
		x->featureInput[i] = 0.0;
    
	for(i=0; i<x->numInstances; i++)
	{
		x->clusterMembers[i].member = (int *)t_getbytes_(2 * sizeof(int));
        
		x->instanceClusterMembership[i]=i; // init membership to index
		x->clusterMembers[i].member[0]=i; // first member of the cluster is the instance index
		x->clusterMembers[i].member[1] = -1;
        
		x->clusterMemberLengths[i] = 2;
	}
    
    
    for(i=0; i<x->numInstances; i++)
    {
		fp = x->instances[i].instance;
        
		for(j=0; j<x->featureLength; j++, fp++)
			fscanf(fd, "%f", fp);
    };
    
    post("read %i instances from file: %s.\n", x->numInstances, buf);    
    post("feature length: %i.", x->featureLength);
	post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
    t_freebytes_(header, 2*sizeof(int));
}


static void timbreID_ARFF(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	FILE *fd;
    int i, j, features_written, att_range_low, att_range_hi;
    t_float *fp;
	t_symbol *filename_symbol, *relation_symbol, *att_symbol;
    char *buf, *filename, *relation, *att_name;
	
	//s=s;
	
	att_range_low = 0;
	att_range_hi = -1;
	att_symbol = 0;
	att_name = 0;
    
	filename_symbol = atom_getsym(argv);
	filename = filename_symbol->s_name;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	if(argc>1)
	{
		relation_symbol = atom_getsym(argv+1);
		relation = relation_symbol->s_name;
	}
	else
		relation = "relation";
    
	fprintf(fd, "@RELATION %s\n\n\n", relation);
    
	if(argc>2)
	{
	    for(i=2; i<argc; i++)
		{
            
			switch((i-2)%3)
			{
				case 0:
					att_range_low = atom_getfloat(argv+i);
					break;
				case 1:
					att_range_hi = atom_getfloat(argv+i);
					break;
				case 2:
					att_symbol = atom_getsym(argv+i);
					att_name = att_symbol->s_name;
					for(j=0; j<=att_range_hi-att_range_low; j++)
						fprintf(fd, "@ATTRIBUTE %s-%i NUMERIC\n", att_name, j);
					break;
				default:
					break;
			}
            
		}
        
		// in case the argument list was incomplete
		for(j=0; j<(x->featureLength-1-att_range_hi); j++)
			fprintf(fd, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", j);
        
	}
	else
	{
		for(i=0; i<x->featureLength; i++)
			fprintf(fd, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", i);
	}
    
    
	fprintf(fd, "\n\n");
	fprintf(fd, "@DATA\n\n");
    
    for(i=0; i<x->numInstances; i++)
    {    
		fp = x->instances[i].instance;
        
		features_written=0; // to keep track of each instances no. of features written.
        
		while(1)
		{
			if(features_written++ == (x->instanceFeatureLengths[i]-1))
			{
				fprintf(fd, "%6.20f", *fp++);
				break;
			}
			else
				fprintf(fd, "%6.20f, ", *fp++);
		};
		
		fprintf(fd, "\n");
   	};
    
    post("wrote %i instances to file: %s.", x->numInstances, buf);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_MATLAB(t_timbreID *x, t_symbol *file_symbol, t_symbol *var_symbol)
{
	FILE *fd;
    int i, features_written;
    t_float *fp;
    char *buf, *filename, *varname;
    
	filename = file_symbol->s_name;
	varname = var_symbol->s_name;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	fprintf(fd, "%% name: %s\n", varname);
	fprintf(fd, "%% type: matrix\n");
	fprintf(fd, "%% rows: %i\n", x->numInstances);
	fprintf(fd, "%% columns: %i\n\n", x->featureLength);
    
    for(i=0; i<x->numInstances; i++)
    {    
		fp = x->instances[i].instance;
        
		features_written=0; // to keep track of each instances no. of features written.
        
		while(1)
		{
			if(features_written++ == (x->instanceFeatureLengths[i]-1))
			{
				fprintf(fd, "%6.20f", *fp++);
				break;
			}
			else
				fprintf(fd, "%6.20f, ", *fp++);
		};
		
		fprintf(fd, "\n");
   	};
    
    post("wrote %i instances to file: %s.", x->numInstances, buf);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_clusters(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *fp;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "wb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	// write a header indicating the number of clusters (x->numClusters)
	fwrite(&x->numClusters, sizeof(int), 1, fd);
    
	fp = x->instanceClusterMembership;
    fwrite(fp, sizeof(int), x->numInstances, fd);
    
	fp = x->clusterMemberLengths;
    fwrite(fp, sizeof(int), x->numClusters, fd);
    
    for(i=0; i<x->numClusters; i++)
    {
		fp = x->clusterMembers[i].member;
		fwrite(fp, sizeof(int), x->clusterMemberLengths[i], fd);
   	};
   	
   	
    post("wrote %i clusters to file: %s.", x->numClusters, buf);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_read_clusters(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *fp, header;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);	
    timbreID_makepath(filename, buf);
    
	// free current x->clusterMembers memory for each list
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
    
	fd = fopen(buf, "rb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }	
    
	// read header indicating number of instruments
	fread(&header, sizeof(int), 1, fd);
    
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), header * sizeof(int));
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), header * sizeof(t_member));
    
	x->numClusters = header;
    
	for(i=0; i<x->numClusters; i++)
		x->clusterMembers[i].member = (int *)t_getbytes_(0);
    
	fp = x->instanceClusterMembership;
    fread(fp, sizeof(int), x->numInstances, fd);
    
	fp = x->clusterMemberLengths;
    fread(fp, sizeof(int), x->numClusters, fd);
    
    for(i=0; i<x->numClusters; i++)
    {
		x->clusterMembers[i].member = (int *)t_resizebytes_(x->clusterMembers[i].member, 0, x->clusterMemberLengths[i] * sizeof(int));
        
		fp = x->clusterMembers[i].member;
		fread(fp, sizeof(int), x->clusterMemberLengths[i], fd);
    };
    
    post("read %i clusters from file: %s.\n", x->numClusters, buf);    
    
	// send the cluster list out the 4th outlet
	timbreID_clusters_list(x);
    
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_clusters_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, j, *fp;
    char *filename = s->s_name;
    char *buf;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		
    
	// write a header indicating number of clusters
	fprintf(fd, "%i\n\n", x->numClusters);
    
    for(i=0; i<x->numClusters; i++)
    {
		fp = x->clusterMembers[i].member;
		
		j=0;
		
		while(x->clusterMembers[i].member[j] != -1)
		{
			fprintf(fd, "%i ", *fp++);
			j++;
		};
        
		fprintf(fd, ";\n");
   	};
   	
    post("wrote %i clusters to file: %s.", x->numClusters, buf);
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_read_clusters_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    int i, j, *fp, header;
    char *filename = s->s_name;
    char *buf, semicolon;
    
	// create local memory
	buf = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
    
    //canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);
    timbreID_makepath(filename, buf);
    
	// free current x->clusterMembers memory for each list
	for(i=0; i<x->numClusters; i++)
		t_freebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i]*sizeof(int));
    
    
	fd = fopen(buf, "r");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }
    
	
	// read header indicating number of instruments
	fscanf(fd, "%i", &header);
    
	x->clusterMemberLengths = (int *)t_resizebytes_(x->clusterMemberLengths, x->numClusters * sizeof(int), header * sizeof(int));
	x->clusterMembers = (t_member *)t_resizebytes_(x->clusterMembers, x->numClusters * sizeof(t_member), header * sizeof(t_member));
	
	x->numClusters = header;
    
    for(i=0; i<x->numClusters; i++)
	{
		// can't be bigger than this, and memory can't be added incrementally during the read. have to get a big block here, then shrink later.
    	x->clusterMemberLengths[i] = x->numInstances;
		x->clusterMembers[i].member = (int *)t_getbytes_(x->clusterMemberLengths[i]*sizeof(int));
	}
	
    for(i=0; i<x->numClusters; i++)
    { 	
		fp = x->clusterMembers[i].member;
		j=0;
		
 		while(fscanf(fd, "%i", fp))
 		{
 			j++;
 			fp++;
 		}
 		
 		fscanf(fd, "%c", &semicolon);
 		
 		// terminate the list with -1
 		*fp = -1;
        
		// shrink off the excess
		x->clusterMembers[i].member = (int *)t_resizebytes_(x->clusterMembers[i].member, x->clusterMemberLengths[i] * sizeof(int), (j+1) * sizeof(int));
		x->clusterMemberLengths[i] = (j+1);
    };
    
    for(i=0; i<x->numClusters; i++)
    {
    	int idx;
    	
    	j=0;
    	
    	while(x->clusterMembers[i].member[j] != -1)
    	{
    		idx = x->clusterMembers[i].member[j];
    		x->instanceClusterMembership[idx] = i;
    		j++;
    	};
    };
    
    post("read %i clusters from file: %s.\n", x->numClusters, buf);    
    
    fclose(fd);
    
    // free memory
    t_freebytes_(buf, MAXPDSTRING*sizeof(char));
}

/* ---------------- utility functions ---------------------- */

static void timbreID_sort_knn_info(int k, int numInstances, int prevMatch, t_knn_info *list)
{
	int i, j, top_i, *top_matches;
	
	top_matches = (int *)t_getbytes_(k * sizeof(int));
    
	for(i=0; i<k; i++)
	{
		t_float max_best;
        
		max_best = FLT_MAX;
		top_i = 0;
        
		for(j=0; j<numInstances; j++)
			if(list[j].dist < max_best)
				if(list[j].idx != prevMatch) // doesn't include previous match - this is good
				{
					max_best = list[j].dist;
					top_i = j;
				};
		
		list[top_i].dist = FLT_MAX;
        
		top_matches[i] = list[top_i].idx;
	}
    
	for(i=0; i<k; i++)
	{
		t_knn_info tmp;
        
		tmp = list[i];
		list[i] = list[top_matches[i]];
		list[top_matches[i]] = tmp;
	}
    
	// free memory
	t_freebytes_(top_matches, k*sizeof(int));
    
	// now, the list passed to the function will have the first k elements in order, 
	// and these elements have the lowest distances in the whole list.
}

static void timbreID_sort_float(int n, t_float *list)
{
	int i, j, flag;
	
	for(i=0; i<n; i++)
	{
		flag = 0;
		
		for(j=0; j<(n-1); j++)
		{
			if(list[j] > list[j+1])
			{
				t_float tmp;
				
				flag = 1;
                
				tmp = list[j+1];
				list[j+1] = list[j];
				list[j] = tmp;
			}
		}
        
		if(flag==0)
			break;
	}
}

static t_float timbreID_mean(int num_rows, int column, t_instance *instances, int normal_flag, t_normData *normData)
{
	int i;
	t_float avg, min, denominator;
	
	avg=0;
	
	for(i=0; i<num_rows; i++)
	{
		if(instances[i].instance[column] == FLT_MAX)
			continue;
		else
		{
			if(normal_flag)
			{
				min = normData[column].min;
				denominator = normData[column].denominator;
                
				avg += (instances[i].instance[column] - min) * denominator;
			}
			else
				avg += instances[i].instance[column];
		}
	}
	
	avg /= num_rows;
	
	return(avg);
}

static t_float timbreID_squared_euclid(t_timbreID *x, t_float *v1, t_float *v2)
{
	int i;
	t_float min, max, norm_denominator, sum, dist;
	
	sum=dist=max=0;
	
	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attributeOrder[i]] < x->normData[x->attributeOrder[i]].min)
				min = v1[x->attributeOrder[i]];
			else
				min = x->normData[x->attributeOrder[i]].min;
            
			if(v1[x->attributeOrder[i]] > x->normData[x->attributeOrder[i]].max)
			{
				max = v1[x->attributeOrder[i]];
                
				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->normData[x->attributeOrder[i]].max;
				norm_denominator = x->normData[x->attributeOrder[i]].denominator;
			}	
			
			dist = (  (v1[x->attributeOrder[i]] - min) * norm_denominator  ) - (  (v2[x->attributeOrder[i]] - min) * norm_denominator  );
			sum += dist*dist*x->weights[x->attributeOrder[i]];
		}
		else
		{
			dist = v1[x->attributeOrder[i]] - v2[x->attributeOrder[i]];
			sum += dist*dist*x->weights[x->attributeOrder[i]];
		}
	}
	
	return(sum);
}

static t_float timbreID_manhattan(t_timbreID *x, t_float *v1, t_float *v2)
{
	int i;
	t_float min, max, norm_denominator, sum, dist;
	
	sum=dist=0;
    
	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attributeOrder[i]] < x->normData[x->attributeOrder[i]].min)
				min = v1[x->attributeOrder[i]];
			else
				min = x->normData[x->attributeOrder[i]].min;
            
			if(v1[x->attributeOrder[i]] > x->normData[x->attributeOrder[i]].max)
			{
				max = v1[x->attributeOrder[i]];
                
				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->normData[x->attributeOrder[i]].max;
				norm_denominator = x->normData[x->attributeOrder[i]].denominator;
			}	
			
			dist = (  (v1[x->attributeOrder[i]] - min) * norm_denominator  ) - (  (v2[x->attributeOrder[i]] - min) * norm_denominator  );
		}
		else
			dist = v1[x->attributeOrder[i]] - v2[x->attributeOrder[i]];
        
		sum += fabs(dist) * x->weights[x->attributeOrder[i]];
	}
	
	return(sum);
}

static t_float timbreID_correlation(t_timbreID *x, t_float *v1, t_float *v2)
{
	int i, j, vecLen, vecLenM1;
	t_float min, max, norm_denominator, vecLenRecip, vecLenM1Recip;
	t_float mean1, mean2, std1, std2, covariance, correlation, *vec1Centered, *vec2Centered;
    
	mean1=mean2=std1=std2=covariance=correlation=0;
	min=max=norm_denominator=1;
    
	vecLen = x->attributehi - x->attributelo + 1;
	vecLenM1 = vecLen-1;
	if(vecLen <= 0) vecLen = 1;
	if(vecLenM1 <=0) vecLenM1 = 1;
	vecLenRecip = 1.0/(t_float)vecLen;
	vecLenM1Recip = 1.0/(t_float)vecLenM1;
    
	vec1Centered = (t_float *)t_getbytes_(vecLen * sizeof(t_float));
	vec2Centered = (t_float *)t_getbytes_(vecLen * sizeof(t_float));
	
	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attributeOrder[i]] < x->normData[x->attributeOrder[i]].min)
				min = v1[x->attributeOrder[i]];
			else
				min = x->normData[x->attributeOrder[i]].min;
            
			if(v1[x->attributeOrder[i]] > x->normData[x->attributeOrder[i]].max)
			{
				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				max = v1[x->attributeOrder[i]];
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->normData[x->attributeOrder[i]].max;
				norm_denominator = x->normData[x->attributeOrder[i]].denominator;
			}	
			
			mean1 += (v1[x->attributeOrder[i]] - min) * norm_denominator;
			mean2 += (v2[x->attributeOrder[i]] - min) * norm_denominator;
		}
		else
		{
			mean1 += v1[x->attributeOrder[i]];
			mean2 += v2[x->attributeOrder[i]];
		}
	};
    
	mean1 *= vecLenRecip;
	mean2 *= vecLenRecip;
    
    
	// min, max, and norm_denominator have already been established
	for(i=x->attributelo, j=0; i<= x->attributehi; i++, j++)
	{	
		if(x->normalize)
		{
			vec1Centered[j] = ((v1[x->attributeOrder[i]] -  min) * norm_denominator) - mean1;
			vec2Centered[j] = ((v2[x->attributeOrder[i]] -  min) * norm_denominator) - mean2;
		}
		else
		{
			vec1Centered[j] = v1[x->attributeOrder[i]] - mean1;
			vec2Centered[j] = v2[x->attributeOrder[i]] - mean2;
		}
	};
    
	for(i=0; i<vecLen; i++)
	{
		std1 += vec1Centered[i]*vec1Centered[i];
		std2 += vec2Centered[i]*vec2Centered[i];
	};
	
	std1 *= vecLenM1Recip;
	std2 *= vecLenM1Recip;
    
	// take sqrt to convert variance to standard deviation
	std1 = sqrt(std1);
	std2 = sqrt(std2);
	
	if(std1==0 || std2==0) std1=std2=0; // don't divide by zero below
    
	for(i=0; i<vecLen; i++)
		covariance += vec1Centered[i]*vec2Centered[i];
    
	// covariance usually averaged via N-1, not N
	covariance *= vecLenM1Recip;
	
	// Pearson correlation coefficient
	correlation = covariance/(std1*std2);
    
	// bash to the 0-2 range, then flip sign so that lower is better. this keeps things consistent with other distance metrics.
	correlation += 1;
	correlation *= -1;
	
	// normally in -1 to 1 range, where higher is more simlar
	// bash it to the 0 to 2 range, where lower is more similar
	//correlation = 2-(correlation+1);
    
	// free memory
	t_freebytes_(vec1Centered, vecLen*sizeof(t_float));
	t_freebytes_(vec2Centered, vecLen*sizeof(t_float));
    
	return(correlation);
}

//static void timbreID_makepath(char *filename, char *buf)
//{
//    char *tmpStr;
//    short defaultPathID, idp = -1;
//    
//	// create local memory
//    tmpStr = (char *)t_getbytes_(MAXPDSTRING * sizeof(char));
//    //buf = strcat(x->defaultPathname, "/");
//    
//    defaultPathID = path_getdefault();
//    
//    idp = path_topotentialname(defaultPathID, filename, tmpStr, 0);
//    if (idp == 0)
//    {
//        path_nameconform(tmpStr, buf, PATH_STYLE_NATIVE, PATH_TYPE_PATH);              
//    }
//    else post("Could not find the file %s ", filename);
//    
//    // free memory    
//    t_freebytes_(tmpStr, MAXPDSTRING*sizeof(char));
//}

static void timbreID_makepath(char *filename, char *buf)
{
    path_nameconform(filename, buf, PATH_STYLE_NATIVE, PATH_TYPE_PATH);             
  
}

/* ---------------- END utility functions ---------------------- */


