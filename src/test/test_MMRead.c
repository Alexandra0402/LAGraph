//------------------------------------------------------------------------------
// LAGraph/src/test/test_MMRead.c:  test cases for LAGraph_MMRead
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include <LAGraph.h>
#include <acutest.h>
#include <graph_zachary_karate.h>

#define OK(method) TEST_CHECK (method == 0)

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

int status ;
GrB_Info info ;
char msg [LAGRAPH_MSG_LEN] ;
// LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL, B = NULL ;
GrB_Type atype = NULL, btype = NULL ;
int version [3] ;
const char *date, *name ;
GrB_Index nrows, ncols, nvals ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// test matrices
//------------------------------------------------------------------------------

// FIXME: this is a relative path, but still fragile
#define DATA_DIR "../data/"
#define TEMP_DIR "/tmp/"

typedef struct
{
    GrB_Index nrows ;
    GrB_Index ncols ;
    GrB_Index nvals ;
    const char *type ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = {
//    nrows ncols  nvals type         name
    {    7,    7,    30, "GrB_BOOL",  "A.mtx" },
    {    7,    7,    12, "GrB_INT32", "cover.mtx" },
    { 1138, 1138,  7450, "GrB_BOOL",  "jagmesh7.mtx" },
    {    8,    8,    18, "GrB_BOOL",  "ldbc-cdlp-directed-example.mtx" },
    {    8,    8,    24, "GrB_BOOL",  "ldbc-cdlp-undirected-example.mtx" },
    {   10,   10,    17, "GrB_BOOL",  "ldbc-directed-example-bool.mtx" },
    {   10,   10,    17, "GrB_FP64",  "ldbc-directed-example.mtx" },
    {   10,   10,    17, "GrB_BOOL",  "ldbc-directed-example-unweighted.mtx" },
    {    9,    9,    24, "GrB_BOOL",  "ldbc-undirected-example-bool.mtx" },
    {    9,    9,    24, "GrB_FP64",  "ldbc-undirected-example.mtx" },
    {    9,    9,    24, "GrB_BOOL",  "ldbc-undirected-example-unweighted.mtx"},
    {   10,   10,    30, "GrB_INT64", "ldbc-wcc-example.mtx" },
    {   14,   14,    46, "GrB_FP64",  "LFAT5.mtx" },
    {    6,    6,     8, "GrB_INT64", "msf1.mtx" },
    {    8,    8,    12, "GrB_INT64", "msf2.mtx" },
    {    5,    5,     7, "GrB_INT64", "msf3.mtx" },
    {    8,    8,    28, "GrB_BOOL",  "sample2.mtx" },
    {    8,    8,    12, "GrB_BOOL",  "sample.mtx" },
    {   64,    1,    64, "GrB_INT64", "sources_7.mtx" },
    { 1000, 1000,  3996, "GrB_FP64",  "olm1000.mtx" },
    { 2003, 2003, 83883, "GrB_FP64",  "bcsstk13.mtx" },
    { 2500, 2500, 12349, "GrB_FP64",  "cryg2500.mtx" },
    {    6,    6,    10, "GrB_INT64", "tree-example.mtx" },
    {   67,   67,   294, "GrB_FP64",  "west0067.mtx" },
    {   27,   51,   102, "GrB_FP64",  "lp_afiro.mtx" },
    {   34,   34,   156, "GrB_BOOL",  "karate.mtx" },
    { 0, 0, 0, "", "" },
    } ;

//------------------------------------------------------------------------------
// typename: return the name of a type
//------------------------------------------------------------------------------

const char *typename (GrB_Type type)
{
    if      (type == GrB_BOOL  ) return ("GrB_BOOL") ;
    else if (type == GrB_INT8  ) return ("GrB_INT8") ;
    else if (type == GrB_INT16 ) return ("GrB_INT16") ;
    else if (type == GrB_INT32 ) return ("GrB_INT32") ;
    else if (type == GrB_INT64 ) return ("GrB_INT64") ;
    else if (type == GrB_UINT8 ) return ("GrB_UINT8") ;
    else if (type == GrB_UINT16) return ("GrB_UINT16") ;
    else if (type == GrB_UINT32) return ("GrB_UINT32") ;
    else if (type == GrB_UINT64) return ("GrB_UINT64") ;
    else if (type == GrB_FP32  ) return ("GrB_FP32") ;
    else if (type == GrB_FP64  ) return ("GrB_FP64") ;
    #if 0
    else if (type == GxB_FC32  ) return ("GxB_FC32") ;
    else if (type == GxB_FC64  ) return ("GxB_FC64") ;
    #endif
    TEST_CHECK (false) ;
    return ("(none)") ;
}

//------------------------------------------------------------------------------
// setup: start a test
//------------------------------------------------------------------------------

void setup (void)
{
    printf ("\nsetup: %s\n", __FILE__) ;
    printf ("data is in [%s]\n", DATA_DIR) ;
    OK (LAGraph_Init (msg)) ;
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, version)) ;
    printf ("%s %d.%d.%d (%s)\n", name,
        version [0], version [1], version [2], date) ;
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    printf ("%s %d.%d.%d (%s)\n", name,
        version [0], version [1], version [2], date) ;
//  OK (LAGraph_Delete (&G, msg)) ;
//  TEST_CHECK (G == NULL) ;
    OK (GrB_free (&A)) ;
    TEST_CHECK (A == NULL) ;
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_MMRead:  read a set of matrices, check their stats, and write them out
//------------------------------------------------------------------------------

void test_MMRead (void)
{
    FILE *f = NULL ;
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, DATA_DIR "%s", aname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        OK (GrB_wait (&A)) ;

        //----------------------------------------------------------------------
        // check its stats
        //----------------------------------------------------------------------

        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;
        OK (GrB_Matrix_nvals (&nvals, A)) ;
        TEST_CHECK (nrows == files [k].nrows) ;
        TEST_CHECK (ncols == files [k].ncols) ;
        TEST_CHECK (nvals == files [k].nvals) ;
        OK (GxB_Matrix_type (&btype, A)) ;
        TEST_CHECK (atype == btype) ;
        OK (GxB_print (A, 2)) ;

        const char *tname = typename (atype) ;
        OK (strcmp (tname, files [k].type)) ;

        //----------------------------------------------------------------------
        // write it to /tmp
        //----------------------------------------------------------------------

        snprintf (filename, LEN, TEMP_DIR "%s", aname) ;
        f = fopen (filename, "w") ;
        OK (LAGraph_MMWrite (A, f, msg)) ;
        OK (fclose (f)) ;

        //----------------------------------------------------------------------
        // load it back in again
        //----------------------------------------------------------------------

        f = fopen (filename, "r") ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        OK (fclose (f)) ;
        OK (GrB_wait (&B)) ;
        OK (GxB_print (B, 2)) ;

        //----------------------------------------------------------------------
        // ensure A and B are the same
        //----------------------------------------------------------------------

        TEST_CHECK (atype == btype) ;

        bool A_and_B_are_identical ;
        OK (LAGraph_IsEqual (&A_and_B_are_identical, A, B, NULL, msg)) ;
        TEST_CHECK (A_and_B_are_identical) ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_karate: read in karate graph from a file and compare it known graph
//-----------------------------------------------------------------------------

void test_karate (void)
{
    setup ( ) ;

    FILE *f = fopen (DATA_DIR "karate.mtx", "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    TEST_CHECK (atype == GrB_BOOL) ;
    OK (fclose (f)) ;
    OK (GxB_print (A, 2)) ;

    OK (GrB_Matrix_new (&B, GrB_BOOL, ZACHARY_NUM_NODES, ZACHARY_NUM_NODES)) ;
    OK (GrB_Matrix_build (B, ZACHARY_I, ZACHARY_J, ZACHARY_V,
        ZACHARY_NUM_EDGES, GrB_LOR)) ;
    OK (GxB_print (B, 2)) ;

    bool A_and_B_are_identical ;
    OK (LAGraph_IsEqual (&A_and_B_are_identical, A, B, NULL, msg)) ;
    TEST_CHECK (A_and_B_are_identical) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "MMRead", test_MMRead },
    { "karate", test_karate },
    { NULL, NULL }
} ;

