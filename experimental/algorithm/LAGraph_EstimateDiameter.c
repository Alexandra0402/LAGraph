//------------------------------------------------------------------------------
// LAGraph_EstimateDiameter: Graph Diameter Estimation
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// Takes in a graph and estimates the diameter 
// and optionally also finds pseudo-peripheral nodes of the graph

// Outputs: 
// Diameter returns the estimated diameter of the graph
// If not set to NULL, peripheral will be a vector with n elements
// index i of peripheral is the estimated diameter if it's a pseudo-peripheral node or nothing if not

// Inputs:
// G is the graph to be analyzed
// maxSrcs limits the number of sources used each cycle
// maxLoops limits the number of times the core loop will run if a stable diameter isn't found
// msg is a buffer for error messages

#define LG_FREE_WORK        \
{                           \
    GrB_free (&ecc) ;       \
    GrB_free (&srcs) ;      \
    GrB_free (&level) ;     \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&peri) ;      \
}

#include "LG_internal.h"

int LAGraph_EstimateDiameter
(
    // outputs:
    GrB_Index    *diameter,
    GrB_Vector    *peripheral,
    // inputs:
    const LAGraph_Graph G,
    GrB_Index    maxSrcs,
    GrB_Index    maxLoops,
    char          *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

#if !LAGRAPH_SUITESPARSE
    LG_ASSERT (false, GrB_NOT_IMPLEMENTED) ;
#else

    LG_CLEAR_MSG ;
    GrB_Vector ecc = NULL ;           // the eccentricity of the nodes
    GrB_Vector peri = NULL ;          // vector to store peripheral node status in
    GrB_Index d = 0 ;              // current diameter
    GrB_Index lastd = 0 ;          // previous diameter
    GrB_Vector srcs = NULL ;        // list of current sources
    GrB_Index nsrcs ;               // number of current sources
    GrB_Matrix level = NULL ;       // matrix for msbfs to put level info in

    bool compute_periphery  = (peripheral != NULL) ;
    if (compute_periphery ) (*peripheral) = NULL ;
    bool compute_diameter  = (diameter != NULL) ;
    LG_ASSERT_MSG (compute_diameter, GrB_NULL_POINTER,
        "Diameter destination must be non-NULL") ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    //--------------------------------------------------------------------------
    // get the problem size and cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    
    GrB_Index n;        // number of nodes in the graph
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    GrB_Type int_type = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;

    //--------------------------------------------------------------------------
    // set up the first maxSrcs random nodes
    //--------------------------------------------------------------------------

    // currently just doing the first maxSrcs, consider different randomization
    // check maxSrcs < n
    // printf("Selecting sources \n");
    if (maxSrcs > n){
        nsrcs = n;
    } else {
        nsrcs = maxSrcs;
    }
    GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;
    for (int64_t i = 0; i < nsrcs; i++){
        GRB_TRY (GrB_Vector_setElement (srcs, i, i)) ;
    }

    //--------------------------------------------------------------------------
    // core loop, run until current and previous diameters match or reach given limit
    //--------------------------------------------------------------------------

    GrB_Monoid max = (n > INT32_MAX) ?
            GrB_MAX_MONOID_INT64 : GrB_MAX_MONOID_INT32 ;
    bool incSrcs = false;
    for (int64_t i = 0; i < maxLoops; i++){
        // printf("Start of main loop \n");
        // save previous diameter
        lastd = d;

        // get new diameter - PUT MSBFS IN TRY
        LAGraph_MultiSourceBFS(&level, NULL, G, srcs, msg);
        GRB_TRY (GrB_Vector_new (&ecc, int_type, n)) ;
        GRB_TRY (GrB_reduce(ecc, NULL, NULL, max, level, GrB_DESC_T0T1)) ;
        GRB_TRY (GrB_reduce(&d, NULL, max, ecc, GrB_NULL)) ;

        // check if done
        if (d == lastd){
            incSrcs = true;
            break;
        }
        // printf("Loop midpoint 1 \n");

        // currently brute forcing by looping through, 
        // continue looking for better solution

        // set up source list for next round
        // get the number of peripheral nodes 
        int64_t nperi = 0;
        for (int64_t j = 0; j < n; j++){
            GrB_Index e;
            GRB_TRY(GrB_Vector_extractElement(&e, ecc, j));
            if (e == d){
                nperi += 1;
            }
        }
        // select the number of sources
        if (nperi > maxSrcs) {
            nsrcs = maxSrcs;
        } else {
            nsrcs = nperi;
        }
        // printf("Loop midpoint 2 \n");
        // printf("Number of peripheral nodes: %d \n",nperi);
        // choose sources
        GrB_free (&srcs) ;
        GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;
        int64_t curri = 0;
        for (int64_t j = 0; j < n; j++){
            GrB_Index e;
            GRB_TRY(GrB_Vector_extractElement(&e, ecc, j));
            if (e == d){
                GRB_TRY (GrB_Vector_setElement (srcs, j, curri)) ;
                curri += 1;
                if (curri == nsrcs) {
                    break;
                }
            }
        }


    }
    // printf("Loop complete \n");

    //--------------------------------------------------------------------------
    // after loop, set up peripheral nodes if needed
    //--------------------------------------------------------------------------
    if (compute_periphery) {
        GRB_TRY (GrB_Vector_new (&peri, int_type, n)) ;
        if (incSrcs) {
            for (int64_t i = 0; i < nsrcs; i++) {
                GrB_Index currsrc;
                GRB_TRY(GrB_Vector_extractElement(&currsrc, srcs, i));
                GRB_TRY (GrB_Vector_setElement (peri, 1, currsrc)) ;  
            }
        }
       
        GRB_TRY (GrB_select(peri, NULL, NULL, GrB_VALUEEQ_T, ecc, d, NULL)) ;

    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    if (compute_periphery) (*peripheral) = peri ;
    (*diameter ) = d ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#endif
}
