//------------------------------------------------------------------------------
// LG_internal.h: include file for use within LAGraph itself
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// These definitions are not meant for the end-user of LAGraph or GraphBLAS

#ifndef LG_INTERNAL_H
#define LG_INTERNAL_H

//------------------------------------------------------------------------------
// include files
//------------------------------------------------------------------------------

#define LG_LIBRARY
#include "LAGraph.h"

#if defined ( __linux__ )
#include <malloc.h>
#endif

#define LG_RESTRICT LAGRAPH_RESTRICT

//------------------------------------------------------------------------------
// string macros
//------------------------------------------------------------------------------

#define LG_XSTR(x) LG_STR(x)
#define LG_STR(x) #x

//------------------------------------------------------------------------------
// string matching
//------------------------------------------------------------------------------

#define MATCH(s1,s2,n) (strncmp (s1, s2, n) == 0)
#define MATCHNAME(s1,s2) MATCH (s1, s2, LAGRAPH_MAX_NAME_LEN)

//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------

// LG_void: used in place of (void *), but valid for pointer arithmetic
typedef unsigned char LG_void ;

//------------------------------------------------------------------------------
// LG_CLEAR_MSG: clear the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method starts, it first clears the caller's msg string.
#define LG_CLEAR_MSG                    \
{                                       \
    if (msg != NULL) msg [0] = '\0' ;   \
}

//------------------------------------------------------------------------------
// LG_ERROR_MSG: set the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method encounters an error, it can report details in the
// msg.  This is normally done via LG_ASSERT_MSG.  For example:

/*
    if (src < 0 || src >= n)
    {
        LG_ERROR_MSG ("Source node %ld must be in range 0 to n-1, "
            "where n = %ld is the number of nodes in the graph.", src, n) ;
        return (GrB_INVALID_INDEX) ;
    }
    // or, with a simpler message:
    LG_ASSERT_MSG (src >= 0 && src < n, GrB_INVALID_INDEX, "invalid src node") ;
*/

#define LG_ERROR_MSG(...)                                           \
{                                                                   \
    if (msg != NULL && msg [0] == '\0')                             \
    {                                                               \
        snprintf (msg, LAGRAPH_MSG_LEN, __VA_ARGS__) ;              \
    }                                                               \
}

//------------------------------------------------------------------------------
// LG_FREE_WORK: free all workspace
//------------------------------------------------------------------------------

#ifndef LG_FREE_WORK
#define LG_FREE_WORK ;
#endif

//------------------------------------------------------------------------------
// LG_FREE_ALL: free all workspace and all output arguments, on error
//------------------------------------------------------------------------------

#ifndef LG_FREE_ALL
#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
}
#endif

//------------------------------------------------------------------------------
// GrB_CATCH: catch an error from GraphBLAS
//------------------------------------------------------------------------------

// A simple GrB_CATCH macro to be used by GrB_TRY.  If an LAGraph function
// wants something else, then #define a GrB_CATCH macro before the #include
// "LG_internal.h" statement.

#ifndef GrB_CATCH

    #define GrB_CATCH(info)                                 \
    {                                                       \
        LG_ERROR_MSG ("%s, line %d: GrB failure: %d",       \
            __FILE__, __LINE__, info) ;                     \
        LG_FREE_ALL ;                                       \
        return (info) ;                                     \
    }

#endif

//------------------------------------------------------------------------------
// LAGraph_CATCH: catch an error from LAGraph
//------------------------------------------------------------------------------

// A simple LAGraph_CATCH macro to be used by LAGraph_TRY.  If an LAGraph
// function wants something else, then #define a LAGraph_CATCH macro before the
// #include "LG_internal.h" statement.

#ifndef LAGraph_CATCH
#define LAGraph_CATCH(status)                           \
{                                                       \
    LG_ERROR_MSG ("%s, line %d: LAGraph failure: %d",   \
        __FILE__, __LINE__, status) ;                   \
    LG_FREE_ALL ;                                       \
    return (status) ;                                   \
}
#endif

//------------------------------------------------------------------------------
// LG_ASSERT_MSG: assert an expression is true, and return if it is false
//------------------------------------------------------------------------------

// Identical to LG_ASSERT, except this allows a different string to be
// included in the message.

#define LG_ASSERT_MSG(expression,error_status,expression_message)       \
{                                                                       \
    if (!(expression))                                                  \
    {                                                                   \
        LG_ERROR_MSG ("LAGraph assertion \"" expression_message         \
            "\" failed:\nfile \"%s\", line %d\n", __FILE__, __LINE__) ; \
        LG_FREE_ALL ;                                              \
        return (error_status) ;                                         \
    }                                                                   \
}

//------------------------------------------------------------------------------
// LG_ASSERT: assert an expression is true, and return if it is false
//------------------------------------------------------------------------------

// LAGraph methods can use this assertion macro for simple errors, such as
// when running of memory:
//
//      void P = LAGraph_Malloc ( ... ) ;
//      LG_ASSERT (P != NULL, GrB_OUT_OF_MEMORY) ;
//
// If LAGraph_Malloc fails and returns P as NULL, the msg is set to:
//
//      LAGraph assertion "P != NULL" failed:
//      file: LAGraph_something, line: 42

#define LG_ASSERT(expression, error_status)                             \
{                                                                       \
    if (!(expression))                                                  \
    {                                                                   \
        LG_ERROR_MSG ("LAGraph assertion \"" LG_XSTR(expression)        \
            "\" failed:\nfile \"%s\", line %d\n", __FILE__, __LINE__) ; \
        LG_FREE_ALL ;                                              \
        return (error_status) ;                                         \
    }                                                                   \
}

//------------------------------------------------------------------------------
// LG_TRY: check a condition and return on error
//------------------------------------------------------------------------------

// The msg is not modified.  This should be used when an LAGraph method calls
// another one.

#define LG_TRY(LAGraph_method)                  \
{                                               \
    int LAGraph_status = LAGraph_method ;       \
    if (LAGraph_status < 0)                     \
    {                                           \
        LG_FREE_ALL ;                      \
        return (LAGraph_status) ;               \
    }                                           \
}

//------------------------------------------------------------------------------
// LG_CLEAR_MSG_AND_BASIC_ASSERT: clear msg and do basic tests of a graph
//------------------------------------------------------------------------------

#define LG_CLEAR_MSG_AND_BASIC_ASSERT(G,msg)                                \
{                                                                           \
    LG_CLEAR_MSG ;                                                          \
    LG_ASSERT (G != NULL, GrB_NULL_POINTER) ;                               \
    LG_ASSERT_MSG (G->A != NULL, LAGRAPH_INVALID_GRAPH,                     \
        "graph adjacency matrix is NULL") ;                                 \
    LG_ASSERT_MSG (G->kind >= LAGRAPH_ADJACENCY_UNDIRECTED &&               \
        G->kind <= LAGRAPH_ADJACENCY_DIRECTED,                              \
        LAGRAPH_INVALID_GRAPH, "graph kind invalid") ;                      \
}

//------------------------------------------------------------------------------
// FPRINTF: fprintf and check result
//------------------------------------------------------------------------------

#define FPRINTF(f,...)                                                      \
{                                                                           \
    LG_ASSERT_MSG (fprintf (f, __VA_ARGS__) >= 0,                           \
        LAGRAPH_IO_ERROR, "Unable to write to file") ;                      \
}

//------------------------------------------------------------------------------
// code development settings
//------------------------------------------------------------------------------

// turn off debugging; do not edit these three lines
#ifndef NDEBUG
#define NDEBUG
#endif

// These flags are used for code development.  Uncomment them as needed.

// to turn on debugging, uncomment this line:
// #undef NDEBUG

#undef ASSERT

#ifndef NDEBUG

    // debugging enabled
    #ifdef MATLAB_MEX_FILE
        // debugging when LAGraph is part of a mexFunction
        #define ASSERT(x)                                               \
        {                                                               \
            if (!(x)) mexErrMsgTxt ("failure: " __FILE__ " line: ") ;   \
        }
    #else
        #include <assert.h>
        #define ASSERT(x) assert (x) ;
    #endif

#else

    // debugging disabled
    #define ASSERT(x)

#endif

//------------------------------------------------------------------------------
// LG_Multiply_size_t:  c = a*b but check for overflow
//------------------------------------------------------------------------------

static bool LG_Multiply_size_t  // true if ok, false if overflow
(
    size_t *c,                  // c = a*b, or zero if overflow occurs
    const size_t a,
    const size_t b
)
{

    ASSERT (c != NULL) ;

    (*c) = 0 ;
    if (a == 0 || b == 0)
    {
        return (true) ;
    }

    if (a > SIZE_MAX / 2 || b > SIZE_MAX / 2)
    {
        // a or b are out of range
        return (false) ;
    }

    // a + b is now safe to compute
    if ((a + b) > (SIZE_MAX / LAGraph_MIN (a,b)))
    {
        // a * b may overflow
        return (false) ;
    }

    // a * b will not overflow
    (*c) = a * b ;
    return (true) ;
}

//------------------------------------------------------------------------------
// Matrix Market format
//------------------------------------------------------------------------------

// %%MatrixMarket matrix <fmt> <type> <storage> uses the following enums:

typedef enum
{
    MM_coordinate,
    MM_array,
}
MM_fmt_enum ;

typedef enum
{
    MM_real,
    MM_integer,
    MM_complex,
    MM_pattern
}
MM_type_enum ;

typedef enum
{
    MM_general,
    MM_symmetric,
    MM_skew_symmetric,
    MM_hermitian
}
MM_storage_enum ;

// maximum length of each line in the Matrix Market file format

// The MatrixMarket format specificies a maximum line length of 1024.
// This is currently sufficient for GraphBLAS but will need to be relaxed
// if this function is extended to handle arbitrary user-defined types.
#define MMLEN 1024
#define MAXLINE MMLEN+6

//------------------------------------------------------------------------------
// LG_PART and LG_PARTITION: definitions for partitioning an index range
//------------------------------------------------------------------------------

// LG_PART and LG_PARTITION:  divide the index range 0:n-1 uniformly
// for nthreads.  LG_PART(tid,n,nthreads) is the first index for thread tid.
#define LG_PART(tid,n,nthreads)  \
    (((tid) * ((double) (n))) / ((double) (nthreads)))

// thread tid will operate on the range k1:(k2-1)
#define LG_PARTITION(k1,k2,n,tid,nthreads)                                  \
    k1 = ((tid) ==  0          ) ?  0  : LG_PART ((tid),  n, nthreads) ;    \
    k2 = ((tid) == (nthreads)-1) ? (n) : LG_PART ((tid)+1,n, nthreads)

//------------------------------------------------------------------------------
// LG_eslice: uniform partition of e items to each task
//------------------------------------------------------------------------------

static inline void LG_eslice
(
    int64_t *Slice,         // array of size ntasks+1
    int64_t e,              // number items to partition amongst the tasks
    const int ntasks        // # of tasks
)
{
    Slice [0] = 0 ;
    for (int tid = 0 ; tid < ntasks ; tid++)
    {
        Slice [tid] = LG_PART (tid, e, ntasks) ;
    }
    Slice [ntasks] = e ;
}

//------------------------------------------------------------------------------
// definitions for sorting functions
//------------------------------------------------------------------------------

// All of the LG_qsort_* functions are single-threaded, by design.  The
// LAGraph_Sort* functions are parallel.  None of these sorting methods are
// guaranteed to be stable.  These functions are contributed by Tim Davis, and
// are derived from SuiteSparse:GraphBLAS v4.0.3.  Functions named LG_* are not
// meant to be accessible by end users of LAGraph.

#define LG_BASECASE (64 * 1024)

void LG_qsort_1a    // sort array A of size 1-by-n
(
    int64_t *LG_RESTRICT A_0,       // size n array
    const int64_t n
) ;

void LG_qsort_2     // sort array A of size 2-by-n, using 2 keys (A [0:1][])
(
    int64_t *LG_RESTRICT A_0,       // size n array
    int64_t *LG_RESTRICT A_1,       // size n array
    const int64_t n
) ;

void LG_qsort_3     // sort array A of size 3-by-n, using 3 keys (A [0:2][])
(
    int64_t *LG_RESTRICT A_0,       // size n array
    int64_t *LG_RESTRICT A_1,       // size n array
    int64_t *LG_RESTRICT A_2,       // size n array
    const int64_t n
) ;

//------------------------------------------------------------------------------
// LG_lt_1: sorting comparator function, one key
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of one integer.

// LG_lt_1 returns true if A [a] < B [b], for LG_qsort_1b

#define LG_lt_1(A_0, a, B_0, b) (A_0 [a] < B_0 [b])

//------------------------------------------------------------------------------
// LG_lt_2: sorting comparator function, two keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of two integers.

// LG_lt_2 returns true if A [a] < B [b], for LG_qsort_2 and LG_msort_2b

#define LG_lt_2(A_0, A_1, a, B_0, B_1, b)                                   \
(                                                                           \
    (A_0 [a] < B_0 [b]) ?                                                   \
    (                                                                       \
        true                                                                \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        (A_0 [a] == B_0 [b]) ?                                              \
        (                                                                   \
            /* primary key is the same; tie-break on the 2nd key */         \
            (A_1 [a] < B_1 [b])                                             \
        )                                                                   \
        :                                                                   \
        (                                                                   \
            false                                                           \
        )                                                                   \
    )                                                                       \
)

//------------------------------------------------------------------------------
// LG_lt_3: sorting comparator function, three keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of three integers.

// LG_lt_3 returns true if A [a] < B [b], for LG_qsort_3 and LG_msort_3b

#define LG_lt_3(A_0, A_1, A_2, a, B_0, B_1, B_2, b)                         \
(                                                                           \
    (A_0 [a] < B_0 [b]) ?                                                   \
    (                                                                       \
        true                                                                \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        (A_0 [a] == B_0 [b]) ?                                              \
        (                                                                   \
            /* primary key is the same; tie-break on the 2nd and 3rd key */ \
            LG_lt_2 (A_1, A_2, a, B_1, B_2, b)                              \
        )                                                                   \
        :                                                                   \
        (                                                                   \
            false                                                           \
        )                                                                   \
    )                                                                       \
)

//------------------------------------------------------------------------------
// LG_eq_*: sorting comparator function, three keys
//------------------------------------------------------------------------------

// A [a] and B [b] are keys of two or three integers.
// LG_eq_* returns true if A [a] == B [b]

#define LG_eq_3(A_0, A_1, A_2, a, B_0, B_1, B_2, b)                         \
(                                                                           \
    (A_0 [a] == B_0 [b]) &&                                                 \
    (A_1 [a] == B_1 [b]) &&                                                 \
    (A_2 [a] == B_2 [b])                                                    \
)

#define LG_eq_2(A_0, A_1, a, B_0, B_1, b)                                   \
(                                                                           \
    (A_0 [a] == B_0 [b]) &&                                                 \
    (A_1 [a] == B_1 [b])                                                    \
)

#define LG_eq_1(A_0, a, B_0, b)                                             \
(                                                                           \
    (A_0 [a] == B_0 [b])                                                    \
)

//------------------------------------------------------------------------------
// count entries on the diagonal of a matrix
//------------------------------------------------------------------------------

int LG_ndiag
(
    // output
    int64_t *ndiag,         // # of entries
    // input
    GrB_Matrix A,           // matrix to count
    char *msg               // error message
) ;

//------------------------------------------------------------------------------
// simple and portable random number generator (internal use only)
//------------------------------------------------------------------------------

#define LG_RANDOM15_MAX 32767
#define LG_RANDOM60_MAX ((1ULL << 60) -1)

// return a random number between 0 and LG_RANDOM15_MAX
GrB_Index LG_Random15 (uint64_t *seed) ;

// return a random uint64_t, in range 0 to LG_RANDOM60_MAX
GrB_Index LG_Random60 (uint64_t *seed) ;

#endif
