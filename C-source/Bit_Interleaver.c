/* file: Bit_Interleaver.c

   Interleaves the input bits according to an input vector that maps the 
   input to the output. The mapping vector is 0-indexed, not 1-indexed, so
   no coversion in C is needed.
   The interleave/deinterleave flag signals whether you want to interleave 
   or deinterleave the input. The possible values are:
     1: interleave
     0: deinterleave

MEXed to LTE_common_bit_interleaver

   The calling syntax is:

      [output] = Bit_Interleaver(input, mapping_vector, interleave/deinterleave )

      output = interleaved input

   Josep Colom Ikuno, INTHFT
   jcolom@nt.tuwien.ac.at
   www.nt.tuwien.ac.at

*/

#include <mex.h>
#include <matrix.h>
#include <stdlib.h>

/* Input Arguments */
#define INPUT       prhs[0]
#define MAPPING     prhs[1]
#define INTERLEAVE  prhs[2]

/* Output Arguments */
#define OUTPUT      plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
    unsigned char   *input_p;
    unsigned char	*output_p;
    double          *mapping;
    mwSize          DataLength;
    mwSize          i;
    mwSize          index;
    double          *interleave;
    
    /* Check if the input is of the correct type */
    if ( (mxIsLogical(INPUT)!=1) && (mxIsClass(INPUT,"uint8")!=1) ) {
        mexErrMsgTxt("Input must be logical or uint8.");
    }
    
    /* Check for proper number of arguments */
    if ((nrhs < 3 )||(nlhs  > 1)) {
        mexErrMsgTxt("Usage: [output] = Bit_Interleaver(input, mapping_vector, interleave/deinterleave)");
    } else {
        /* first input is the data word */
        input_p = mxGetPr(INPUT);
        mapping = mxGetPr(MAPPING);
        interleave = mxGetPr(INTERLEAVE);
        DataLength = mxGetN(INPUT); /* number of data bits */
        
        /* create the output vector */
        if(mxIsLogical(INPUT)==1) {
            OUTPUT = mxCreateLogicalMatrix(1, DataLength);
        } else {
            OUTPUT = mxCreateNumericMatrix(1, DataLength,mxUINT8_CLASS,mxREAL);
        }
        output_p = mxGetPr(OUTPUT);
        
        /* Populate the output */
        for(i=0;i<DataLength;i++) {
            index = (int)mapping[i];
            if(*interleave==1) output_p[i] = input_p[index];
            else  output_p[index] = input_p[i];
        }
    }
    
    return;
}
