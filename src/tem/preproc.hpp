#ifndef PREPROC_H
#define PREPROC_H

// When using the Borland compiler DEFINE the compiler directive 
//   BORLAND_CPP below. Otherwise, DEFINE the compiler 
//   directive ANSI_CPP below. 

// NOTE:  The Borland C++ Builder compiler stores time variables
//        in a different library <time> than the g++ or
//        Portland Group C++ complier <ctime>.  When using the 
//        Borland compiler DEFINE the compiler directive 
//        BORLAND_CPP below. Otherwise, DEFINE the compiler 
//        directive ANSI_CPP below  

//#define BORLAND_CPP
#define ANSI_CPP


// NOTE: If running TEM "offline", make sure the compiler 
//         directive STANDALONE_TEM is DEFINED below.  If  
//         coupling TEM to the IGSM, make sure STANDALONE_TEM is 
//         NOT DEFINED below (i.e. commented out)

//#define STANDALONE_TEM

#endif
