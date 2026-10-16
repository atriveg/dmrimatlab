/*==========================================================
 * carlsonR_.cpp
 *
 * Used to evaluate elliptic integrals of several kinds
 * over vectors x, y, and z (parallel implementation)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2026 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/CarlsonR.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"


typedef struct CarlsonRIOData{
    BufferType x;
    BufferType y;
    BufferType z;
    BufferType s;
    BufferType r;
} CarlsonRIOData;

typedef struct CarlsonParameters{
    unsigned char kind;
    unsigned int  singletons[3];
} CarlsonParameters;

class ThArgs : public DMRIThreader
{
public:
    CarlsonRIOData*    io;
    CarlsonParameters* params;
};

THFCNRET carlsonR_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     *
     * prhs[0]:  x, the vector x, Nx1
     * prhs[1]:  y, the vector y, Nx1
     * prhs[2]:  z, the vector z, Nx1
     * prhs[3]:  kind, either 'd' or 'f', the kind of elliptic integral to compute
     * prhs[4]:  maxthreads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: the value of Rd(x,y,z) or Rf(x,y,z)
     * plhs[1]: status information about the computation
     *
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:carlsonR_:callstack","This function should only be called from carlsonRd or carlsonRf");
    else if( strcmp(callerFunc,"carlsonRd") && strcmp(callerFunc,"carlsonRf") )
        mexErrMsgIdAndTxt("MyToolbox:carlsonR_:callstack","This function should only be called from carlsonRd or carlsonRf");
    //=======================================================================================
    if(nrhs!=5)
        mexErrMsgIdAndTxt("MyToolbox:carlsonR_:nrhs","Exactly 3 input arguments are required");
    if(nlhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:carlsonR_:nlhs","Exactly 2 output argument are required");
    //=======================================================================================
    // ------ Inputs
    CarlsonRIOData io;
    io.x = mxGetDoubles(prhs[0]);
    io.y = mxGetDoubles(prhs[1]);
    io.z = mxGetDoubles(prhs[2]);
    // ------ Parameters
    CarlsonParameters params;
    params.kind = mxGetChars(prhs[3])[0];
    SizeType N = 1;
    for( unsigned int d=0; d<3; ++d ){
        SizeType Nd = mxGetM(prhs[d]);
        if( Nd==1 )
            params.singletons[d] = 0;
        else{
            params.singletons[d] = 1;
            N = Nd;
        }
    }
    // ------ Outputs
    plhs[0] = mxCreateDoubleMatrix( N, 1, mxREAL );
    io.s = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( N, 1, mxREAL );
    io.r = mxGetDoubles(plhs[1]);
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[4]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.io         = &io;
    threader.params     = &params;
    //=======================================================================================
    threader.threadedProcess( maxthreads, carlsonR_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET carlsonR_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    CarlsonRIOData* io = args->io;
    CarlsonParameters* params = args->params;

    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for(  IndexType i=start; i<end; ++i ){
            if( params->kind=='d'){
                io->r[i] = ElementType(   carlson::CarlsonRd(
                    io->x[i*(params->singletons[0])],
                    io->y[i*(params->singletons[1])],
                    io->z[i*(params->singletons[2])], io->s[i]
                                                            )
                                      );
            }
            else{ // i.e. params->kind == 'f'
                io->r[i] = ElementType(   carlson::CarlsonRf(
                    io->x[i*(params->singletons[0])],
                    io->y[i*(params->singletons[1])],
                    io->z[i*(params->singletons[2])], io->s[i]
                                                            )
                                      );
            }
        }
    }
    while( start < args->getN() );

    return (THFCNRET)NULL;
}
