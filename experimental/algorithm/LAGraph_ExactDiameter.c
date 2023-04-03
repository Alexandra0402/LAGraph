//------------------------------------------------------------------------------
// LAGraph_ExactDiameter: Graph Diameter Computation
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// Takes in a graph and finds the diameter 
// and optionally also the peripheral nodes of the graph

#define LG_FREE_WORK        \
{                           \
    GrB_free (&ecc) ;       \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&peri) ;      \
}

#include "LG_internal.h"

int LAGraph_ExactDiameter
(
    // outputs:
    GrB_Scalar    *diameter,
    GrB_Vector    *peripheral,
    // inputs:
    const LAGraph_Graph G,
    GrB_Index      k,
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
    GrB_Index d ;                     // diameter

    bool compute_periphery  = (peripheral != NULL) ;
    if (compute_periphery ) (*peripheral) = NULL ;
    bool compute_diameter  = (diameter != NULL) ;
    if (compute_diameter ) (*diameter) = NULL ;
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
    GRB_TRY (GrB_Vector_new (&ecc, int_type, n)) ;

    //--------------------------------------------------------------------------
    // get eccentricity, k nodes at a time
    //--------------------------------------------------------------------------

    int64_t setStart = 0;
    GrB_Monoid max = (n > INT32_MAX) ?
            GrB_MAX_MONOID_INT64 : GrB_MAX_MONOID_INT32 ;
    while (setStart < n){
        // set up the sources for this iteration
        GrB_Vector srcs;
        int64_t nsrcs;
        if ((setStart + k) <= n){
            nsrcs = k;
        } else {
            nsrcs = n - setStart;
        }
        GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;
        GrB_Index sources[nsrcs];
        for (int64_t i = 0; i < nsrcs; i++){
            GRB_TRY (GrB_Vector_setElement (srcs, setStart+i, i)) ;
            sources[i] = setStart+i;
        }

        // run bfs to get level matrix for the sources
        GrB_Matrix level = NULL;
        LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, NULL, G, srcs, msg)) ;

        // populate setStart to setStart+nsrcs of ecc with the max level for each src
        GrB_Vector srcEcc;
        GRB_TRY (GrB_Vector_new (&srcEcc, int_type, nsrcs)) ;
        GRB_TRY (GrB_reduce(srcEcc, NULL, NULL, max, level, GrB_NULL)) ;
        LAGRAPH_TRY (GxB_subassign(ecc, NULL, NULL, srcEcc, sources,  nsrcs, GrB_NULL)) ;

        // adjust setStart for next iteration
        setStart = setStart + nsrcs;
    }


    //--------------------------------------------------------------------------
    // determine diameter from eccentricity list
    //--------------------------------------------------------------------------
    
    GRB_TRY (GrB_reduce(&d, NULL, max, ecc, GrB_NULL)) ;

    //--------------------------------------------------------------------------
    // get peripheral nodes, if applicable
    //--------------------------------------------------------------------------

    // currently brute forcing by looping through, 
    // continue looking for better solution 
    if (compute_periphery){
        GRB_TRY (GrB_Vector_new (&peri, int_type, n)) ;
        for (int64_t i = 0; i < n; i++){
            GrB_Index e;
            GRB_TRY(GrB_Vector_extractElement(&e, ecc, i));
            if (e == d){
                GRB_TRY (GrB_Vector_setElement (peri, 1, i)) ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    if (compute_periphery) (*peripheral) = peri ;
    GRB_TRY(GrB_Scalar_new(diameter, int_type)) ;
    GRB_TRY(GrB_Scalar_setElement(*diameter, d)) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#endif
}
