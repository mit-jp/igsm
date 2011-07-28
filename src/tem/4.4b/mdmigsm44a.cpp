/* *************************************************************
****************************************************************
MDMIGSM44A.CPP - Methane Dynamics Model 

Modifications
20021231 - Daily version created by Q. Zhuang
20030103 - Q. Zhuang added new hydrological model
20080130 - DWK changed include from mdm51.h to mdmigsm44a.h
20080130 - DWK changed MDM51:: to MDM44::
 
****************************************************************
************************************************************* */

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iostream>

  using std::cout;
  using std::endl;

#include<cstdlib>

  using std::exit;

#include<cmath>

  using std::ceil;

#include<vector>

  using std::vector;
    
#include<string>
  
  using std::string;


#include "mdmigsm44a.h"


/* ********************************************************** */

MDM44::MDM44() : predstr( NUMMDM )
{
  // Identify potential output variables from MDM
  
  // Gross amount of atmospheric methane taken up by soil 
  predstr[I_CH4CON] = "CH4CON";

  // Gross amount of methane emitted from soil
  predstr[I_CH4EMI] = "CH4EMI";

  // Net amount of methane emitted (or consumed) from soil 
  predstr[I_TOTFLX] = "CH4TOT";

};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void MDM44::stepday( const int& pcmnt,
                     const int& picult,
                     const int& pwet,
                     const double& pfrin,
                     const double& ph2otable,
                     const double& ppctsand,
                     const double& ppctsilt,
                     const double& ppctclay,
                     const double& ppcfldcap,
                     const double& prootz,
                     const double& pph,
                     const double& psuborg,
                     const double& ptsoil )
{
  int i;
  
  double depthz;
  
  int satzone;
    
  double z;

  
  // Initialize CH4 fluxes for each 1 cm depth
  
  for( i = 0; i < MXMDMNLAY; ++i )
  {    
    soil.setCH4RATESAT( ZERO, i );
    soil.setCH4OXIRATE1( ZERO, i );
    soil.setFDFS( ZERO, i );
  }


  // Initialize methane consumption to zero for time step
  
  dryMethanotroph.setCONSP( ZERO );
  wetMethanotroph.setCONSP( ZERO );
  soil.setCH4CONS( ZERO );

  // If entire area is not inundated (i.e. pfrin == 10000.0),
  //   determine methane consumption in upland areas

  if( pfrin < 10000.0 )  // upland
  {
    // Determine methane diffusion into the soil 
    //   (i.e. initially assume no methanotrophy )
    
    soil.setSAT( 0 ); // unsaturated soils
    
    soil.setFDF( soil.setDiffusivityFactor( ppctsand, 
                                            ppctsilt, 
                                            ppctclay, 
                                            soil.getSAT() ) );

    soil.setCH4FLX( soil.dryCH4diffusion( soil.getFDF(), 
                                          methanogen.getLOWB( pcmnt ) ) );
      
    // Determine methanotrophy between surface and low boundary 
      
    for( i = 0; 
         i < (int) (methanogen.getLOWB( pcmnt ) /10.0); 
         ++i )
    {

      soil.setOxygenConc( dryMethanotroph.getAFP( pcmnt ), 
                          methanogen.getLOWB( pcmnt ) );

      depthz = (i+1) * 10.0;
      
      soil.setEHLT( soil.RedoxP( ph2otable, 
                                 depthz, 
                                 ppcfldcap, 
                                 pcmnt, 
                                 methanogen.getLOWB( pcmnt ) ) );


      // Set CH4 oxidation rate equal to the CH4 consumption 
      //   rate by methanotrophs for each 1 cm soil layer
      
      soil.setCH4OXIRATE1( dryMethanotroph.consumeCH4( pcmnt, 
                                                       soil.getCH4CON( i ), 
                                                       soil.getOXYGENC( i ), 
                                                       soil.getINTERSOILT( i ),
                                                       soil.getINTERSH2O( i ),
                                                       soil.getEHLT(), 
                                                       picult ),
                           i );

     
      // Aggregate soil layer CH4 consumption to determine 
      //   total consumption within the whole soil profile 
        
      dryMethanotroph.setCONSP( (dryMethanotroph.getCONSP() 
                                + soil.getCH4OXIRATE1( i )) );
    }
    

    // Determine CH4 uptake into soils (mg m-2 day-1) as a 
    //   result of both diffusion and methanotrophy 
    
    soil.setCH4FLX( soil.dryCH4uptake( soil.getFDF(), 
                                       methanogen.getLOWB( pcmnt ) ) );


    // Restrict CH4 uptake to be greater than or equal to zero
        
    if( soil.getCH4FLX() < ZERO ) { soil.setCH4FLX( ZERO ); }
    
    
    // Assume no CH4 uptake occurs when temperature is less than
    //   -1.0 degree Celsius

    if( ptsoil < -1.0 ) { soil.setCH4FLX( ZERO ); } 


    // Weight upland CH4 uptake by proportion of grid cell
    //   considered to be upland

    soil.setCH4CONS( (soil.getCH4FLX() * (1.0 - (pfrin / 10000.0))) );   

    if( (soil.getCH4CONS() < MISSING) 
        || (soil.getCH4CONS() > 999999.9) )
    {
      soil.setCH4CONS( MISSING );
    }
  }  // end of upland consumption module


  // Calculate methane production rate at wetland 
  //   --- with daily time step

  if( 0 == pwet ) { soil.setCH4EMIS( ZERO ); }
  else // wetland
  {     
    if( ph2otable < ZERO )
    {
      // water table above soil surface, methanogenesis occured 
      //   above the low boundary, oxidation only occured as 
      //   rhizopherix zone 
      // transport will include diffusion and ebullision and 
      //   plant aieded emissions,
    }
    else
    {
      // water table below soil surface
      
      // Determine soil layer containing the water table
      
      satzone = (int) (ceil( ph2otable / 10.0 ));
      
 
      // Determine methanogenesis in saturated soils below the 
      //   water table and above the lower boundary

      for( i = satzone; 
           i < ceil( (methanogen.getLOWB( pcmnt ) / 10.0) ); 
           ++i )
      {               
        depthz = i * 10.0;
        
        soil.setEHLT( soil.RedoxP( ph2otable, 
                                   depthz, 
                                   ppcfldcap, 
                                   pcmnt, 
                                   methanogen.getLOWB( pcmnt ) ) );
       

        // Determine CH4 production within each 1 cm soil layer
        
        soil.setCH4RATE( methanogen.produceCH4( pcmnt, 
                                                psuborg, 
                                                depthz,
                                                prootz,
                                                soil.getINTERSOILT( i ),
                                                pph,
                                                soil.getEHLT() ),
                         i );       
      }

      
      // Determine loss of CH4 from the soil profile due to 
      //   ebullition (i.e. bubble formation) in the 
      //   saturated soil layers between the water table and 
      //   the lower boundary      
     
      soil.setEBUFLX( ZERO );
      
//      for( i = satzone; 
//           i < ceil( (methanogen.getLOWB( pcmnt ) / 10.0) ); 
//           ++i )
//      {

        // Determine CH4 ebullition from each 1 cm soil layer
        
//        soil.setCH4EBUFLX( soil.CH4ebullition( soil.getCH4RATE( i ) ),
//                           i );
     
        // Aggregate soil layer CH4 ebullitions to determine 
        //   total ebullition flux from the whole soil profile 
        //   and covert units from umol L-1 hr-1 to 
        //   mg CH4 m-2 day-1:
        //   CH4EBUFLX = CH4EBUFLX * 0.001 // umol cm-3 hr-1
        //               * 10000.0      // u mol m-2 hr-1
        //               * 0.000001     // mol m-2 hr-1
        //               * 16.0         // g m-2 hr-1
        //               * 24.0         // g m-2 day-1
        //               * 1000.0       // mg m-2 day-1
        //             = CH4EBUFLX * 3.84

//        soil.setEBUFLX( (soil.getEBUFLX() 
//                        + (soil.getCH4EBUFLX( i ) * 3.84)) );
//      }

      
      // Determine loss of CH4 from soil profile due to 
      //   plant-aided transport and rhizospheric CH4 oxidation
    
      soil.setTVEG( 1.0 );
      
      soil.setCH4PLTFLX( ZERO );
   
      for( i = satzone; 
           i < ceil( (methanogen.getLOWB( pcmnt ) / 10.0) ); 
           ++i )
      {
        z = i * 10.0;
      
      
        // Determine the amount of CH4 available for plant 
        //   transport in each 1 cm soil layer by subtracting 
        //   the ebullition rate from the production rate in 
        //   each layer
        
//        soil.setCH4RATE( (soil.getCH4RATE( i )
//                          - soil.getCH4EBUFLX( i )),
//                         i ); 
      
           
        // Determine the potential loss of CH4 from each soil layer
        //   due to plant transport
        
        soil.setCH4PLANTR( soil.CH4plantTransport( soil.getTVEG(), 
                                                   prootz,
                                                   z, 
                                                   soil.getINTERSOILT( i ), 
                                                   soil.getCH4RATE( i ) ),
                           i );

      
        // Update the change of CH4 concentration in each soil 
        //   layer, by subtracting the loss of CH4 due to
        //   plant transport from soil.ch4rate
      
        soil.setCH4RATE( (soil.getCH4RATE( i )
                          - soil.getCH4PLANTR( i )),
                         i );


        // Determine net increase in CH4 concentrations in
        //   each saturated soil layer
         
        soil.setCH4RATESAT( soil.getCH4RATE( i ), i );


        // Assume rhizopheric oxidation of CH4 account for 
        //   40% of ch4plantr[i], the remaining CH4 (60%) is 
        //   transported to the atmosphere   
      
        soil.setCH4PLANTR( (soil.getCH4PLANTR( i ) * 0.60), i );
        

        // Aggregate CH4 flux by plant transport from each soil 
        //   layer to determine total plant transport flux from 
        //   the entire soil profile 
        //   and covert units from umol L-1 hr-1 to 
        //   mg CH4 m-2 day-1:
        //   CH4PLANTR = CH4PLANTR * 0.001 // umol cm-3 hr-1
        //               * 10000.0         // u mol m-2 hr-1
        //               * 0.000001        // mol m-2 hr-1
        //               * 16.0            // g m-2 hr-1
        //               * 24.0            // g m-2 day-1
        //               * 1000.0          // mg m-2 day-1
        //             = CH4PLANTR * 3.84

        soil.setCH4PLTFLX( (soil.getCH4PLTFLX() 
                           + (soil.getCH4PLANTR( i ) * 3.84)) );
      }


      // Determine diffusivity factors for both the saturated 
      //   and unsaturated soil layers in wetlands
            
      for( i = satzone; 
           i < ceil( (methanogen.getLOWB( pcmnt ) / 10.0) ); 
           ++i )
      {
        soil.setSAT( 1 );  // saturated soils
        
        soil.setFDFS( soil.setDiffusivityFactor( ppctsand, 
                                                 ppctsilt, 
                                                 ppctclay, 
                                                 soil.getSAT() ),
                      i );
      }
      
      for( i = 0; i < satzone; ++i )
      {
        soil.setSAT( 0 );  // unsaturated soils
        
        soil.setFDFS( soil.setDiffusivityFactor( ppctsand, 
                                                 ppctsilt, 
                                                 ppctclay, 
                                                 soil.getSAT() ),
                      i );
      }

      // Determine redistribution of CH4 within the soil profile due
      //   to diffusion

      soil.setCH4FLX( soil.wetCH4diffusion( methanogen.getLOWB( pcmnt ),
                                            ph2otable ) );

           
      // Determine soil moisture in unsaturated zone above the
      //   water table (i.e. unsatthetaWL[i]) based on water
      //   table depth
      
      soil.setUnsaturatedZoneMoisture( ph2otable );      

    // Determine methanotrophy in the unsaturated soil layers
    //   between the soil surface and the water table 
      
      for( i = 0; i < satzone; ++i )
      {
        soil.setOxygenConc( wetMethanotroph.getAFP( pcmnt ), 
                            methanogen.getLOWB( pcmnt ) );

        depthz = (i+1) * 10.0;
        
        soil.setEHLT( soil.RedoxP( ph2otable, 
                                   depthz, 
                                   ppcfldcap, 
                                   pcmnt,
                                   methanogen.getLOWB( pcmnt ) ) );


        // Set CH4 oxidation rate within each 1 cm unsaturated 
        //   soil layer equal to the CH4 consumption rate by 
        //   methanotrophs in that layer
       
        soil.setCH4OXIRATE1( wetMethanotroph.consumeCH4( pcmnt,
                                                         soil.getCH4CON_B( i ), 
                                                         soil.getOXYGENC( i ), 
                                                         soil.getINTERSOILT( i ),
                                                         soil.getUNSATTHETAWL( i ),
                                                         soil.getEHLT(), 
                                                         picult ),
                             i );
                                                                    
      }
     
     
      // Determine net CH4 emissions due to diffusion in wetlands
      
      soil.setCH4FLX( soil.wetCH4flux( methanogen.getLOWB( pcmnt ),
                                       ph2otable, 
                                       ptsoil ) );


      // Determine total CH4 net emissions from wetlands from 
      //   diffusion, plant transport and ebullition. 
      //   Weight CH4 emissions by proportion of grid cell
      //   considered to be wetland 

//      soil.setCH4EMIS( ((soil.getCH4FLX() 
//                       - soil.getCH4PLTFLX() 
//                       - soil.getEBUFLX())
//                       * (pfrin / 10000.0)) );

      // Note: consider only diffusion and plant transport

      soil.setCH4EMIS( ((soil.getCH4FLX() 
                       - soil.getCH4PLTFLX())
                       * (pfrin / 10000.0)) );
           
    } // end of else for water table
      

    if( (soil.getCH4EMIS() < MISSING) 
        || (soil.getCH4EMIS() > 999999.9) )
    {
      soil.setCH4EMIS( MISSING );
    }
  } //end for wetland


  // Determine overall CH4 flux for grid cell 
  //   (i.e. CH4 emissions from wetlands plus 
  //   CH4 uptake in uplands - note: CH4 soil emissions occur
  //   with negative values of CH4TOT)

  soil.setCH4TOT( (soil.getCH4EMIS() + soil.getCH4CONS()) );

  if( (soil.getCH4TOT() < MISSING) 
      || (soil.getCH4TOT() > 999999.9) )
  {
    soil.setCH4TOT( MISSING );
  }
	
};



