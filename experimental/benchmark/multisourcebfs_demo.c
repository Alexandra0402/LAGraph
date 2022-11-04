//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/multisourcebfs_demo.c: a simple demo
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGraph_MultiSourceBFS algorithm, in experimental/algorithmn based on 
// helloworld2_demo.  To use it, compile LAGraph while in the build folder
// with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/multisourcebfs_demo < ../data/west0067.mtx
//      ./experimental/benchmark/multisourcebfs_demo < ../data/karate.mtx
//
#include "LAGraphX.h"

// LAGRAPH_CATCH is required by LAGRAPH_TRY.  If an error occurs, this macro
// catches it and takes corrective action, then terminates this program.
#define LAGRAPH_CATCH(info)                     \
{                                               \
    GrB_free (&Y) ;                             \
    GrB_free (&A) ;                             \
    GrB_free (&parent) ;                        \
    GrB_free (&level) ;                         \
    GrB_free (&SourceNodes) ;                   \
    LAGraph_Delete (&G, msg) ;                  \
    return (info) ;                             \
}

// GRB_CATCH is required by GRB_TRY (although GRB_TRY isn't used here)
#define GRB_CATCH(info) LAGRAPH_CATCH(info)

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix Y = NULL, A = NULL ;
    GrB_Matrix level = NULL ;
    GrB_Matrix parent = NULL ;
    GrB_Vector SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    LAGRAPH_TRY (LAGraph_Init (msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph via a Matrix Market file from stdin
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MMRead (&A, stdin, msg)) ;
    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the LAGraph_MultiSourceBFS algorithm
    //--------------------------------------------------------------------------

    printf ("\n==========================Set up for BFS\n") ;
    // set up
    GRB_TRY (GrB_Vector_new (&SourceNodes, GrB_INT32, 1)) ;
    printf ("\n==========================Intermediate Print\n") ;
    GRB_TRY (GrB_Vector_setElement (SourceNodes, 0, 0)) ; // currently just uses node 0 of the graph

    printf ("\n==========================Running BFS\n") ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, &parent, G, SourceNodes, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_MultiSourceBFS: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure Y is a copy of G->A)
    //--------------------------------------------------------------------------
/*
    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_Matrix_IsEqual (&isequal, Y, G->A, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;
    if (isequal)
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }
*/
    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix level:\n") ;
    LAGRAPH_TRY (LAGraph_Matrix_Print (level, LAGraph_SHORT, stdout, msg)) ;
    printf ("\n===============================The result matrix parent:\n") ;
    LAGRAPH_TRY (LAGraph_Matrix_Print (parent, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    GrB_free (&Y) ;
    LAGraph_Delete (&G, msg) ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

