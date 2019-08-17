 /*
 Gold_sequence_generation - generates the gold sequence needed to scramble and descramble the coded bits in each of the codewords for PDSCH.
 DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03) Section 7.2
 gold_seq = Gold_sequence_generation(output_seq_length, Cinit)
 Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at / 2010 modified by Petr Kejik, xkejik00@stud.feec.vutbr.cz
 (c) 2016 by ITC 
 www.nt.tuwien.ac.at

 MEXed to LTE_common_gen_gold_sequence_CCH

 input :   output_seq_length   ... [1 x 1]double - length of coded bits from one codeword = length of gold sequence
 output:   gold_seq            ... [1 x # coded bits]double - the gold sequence for one codeword

 date of creation: 2009/03/19
*/

#include <stdlib.h>
#include <math.h>
#include <mex.h>

/* Input Arguments */
#define INPUT_length     prhs[0]
#define INPUT_Cinit      prhs[1]

/* Output Arguments */
#define OUTPUT_gold_seq  plhs[0]

/* main function that interfaces with MATLAB */
void mexFunction(
				 int            nlhs,
				 mxArray       *plhs[],
				 int            nrhs,
				 const mxArray *prhs[] )
{
  /* polynomials defining two m-sequences */
  /* DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03) Section 7.2 */
  unsigned char poly0[]       = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
  unsigned char poly1[]       = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1};
  
  /* the mask to get the next state of the first m-seqence */
  unsigned char mseq_mask0[]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  /* the mask to get the next state of the second m-seqence */
  unsigned char mseq_mask1[]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  
  /*initial states of m-sequence generators
    DEFINED IN STANDARD 3GPP TS 36.211 V8.2.0 (2008-03) Section 7.2
    PN register state of the first m-seqence */
  unsigned char mseq_state0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  /* PN register state of the second m-seqence */
  unsigned char mseq_state1[31];
  long seed1_dec; /* initial state expressed as decimal value */
  
  mwSize output_seq_length;
  double Cinit;
  unsigned char *gold_seq;

  unsigned char mseq0; /* output of the first m-sequence */
  unsigned char mseq1; /* output of the second m-sequence */
  int deg; /* order of the polynomials poly0 and poly1 */
  unsigned char fdbkBit0; /* feedback bit of the first m-sequence */
  unsigned char fdbkBit1; /* feedback bit of the second m-sequence */
  unsigned char tmp0; /* temporary output variable of the first m-sequence */
  unsigned char tmp1; /* temporary output variable of the second m-sequence */
  
  mwSize i, ji; /* loop variables */
  int help_convert; /* temporary variable for the dec to binary conversion */
  
  /* Check for proper number of arguments */
  if ((nrhs != 2 )||(nlhs  < 1)||(nlhs  > 2)) {
	mexErrMsgTxt("Usage: s = LTE_common_scrambling(b, NIDcell, NIDue, subframe, codeword, mode)");
  } 
  else {
	/* Get input parameters */
	output_seq_length = (mwSize) (*mxGetPr(INPUT_length));
	Cinit             = *mxGetPr(INPUT_Cinit);
    
    /* create the output vector */
    OUTPUT_gold_seq = mxCreateLogicalMatrix(1,output_seq_length);
    gold_seq = mxGetPr(OUTPUT_gold_seq);

    deg = 31;    
    seed1_dec = Cinit;
    
    /* Conversion of the seed into binary number for the second m-sequence that is calculated as decimal from input pars */
    for(help_convert = 30; help_convert >= 0; help_convert--) {
        mseq_state1[help_convert] = seed1_dec / (1 << help_convert);
        seed1_dec = seed1_dec - mseq_state1[help_convert] * (1 << help_convert);
    }
    
    for (i = 0; i < output_seq_length; i++) {
        /* Compute feedback bit */
        fdbkBit0 = 0;
        fdbkBit1 = 0;
        for (ji = 1; ji < deg + 1; ji++) {
            fdbkBit0 = fdbkBit0 ^ (poly0[ji] * mseq_state0[ji-1]);
            fdbkBit1 = fdbkBit1 ^ (poly1[ji] * mseq_state1[ji-1]);
        }
      
        /* Apply output mask and compute output bit */
        tmp0 = 0;
        tmp1 = 0;
        for (ji = 0; ji < deg; ji++) {
            tmp0 = tmp0 ^ (mseq_state0[ji] & mseq_mask0[ji]);
            tmp1 = tmp1 ^ (mseq_state1[ji] & mseq_mask1[ji]);
        }
        mseq0 = tmp0;
        mseq1 = tmp1;
        
        /* Update states */
        for (ji = deg-1; ji >= 1; ji--) {
            mseq_state0[ji] = mseq_state0[ji-1];
            mseq_state1[ji] = mseq_state1[ji-1];
        }
        mseq_state0[0] = fdbkBit0;
        mseq_state1[0] = fdbkBit1;
        gold_seq[i] = mseq0 ^ mseq1;
    }
  }
  return;
}
