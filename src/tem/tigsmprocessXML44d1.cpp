/* *************************************************************
TIGSMPROCESSMXL44C.CPP - Functions to process XML files for TEM

Programmer: David Kicklighter
Creation Date: 17 June 2003

Modifications:

20060126 - DWK created by modifying tprocessXML51.cpp
20060126 - DWK changed include from tprocessXML51.h to 
           tprocessXML437.h
20060126 - DWK changed ProcessXML50:: to ProcessXML43::
20080130 - DWK changed include from tprocessXML437.h to
           tigsmprocessXML44a.h
20080130 - DWK changed ProcessXML43:: to ProcessXML44::
20110707 - DWK changed include from tigsmprocessXML44a.h to
           tigsmprocessXML44c.h
20150429 - DWK changed include tigsmprocessXML44c.h to
           tigsmprocessXML44d1.h
           
************************************************************* */

#include<iostream>

  using std::cerr;
  using std::endl;

#include<fstream>

  using std::ifstream;

#include<cstdlib>

  using std::exit;
  using std::atof;
  using std::atoi;
  using std::atol;
  
#include<string>
  
  using std::string;

#include "tigsmprocessXML44d1.h"


ProcessXML44::ProcessXML44()
{

};

/* *************************************************************
************************************************************* */

void ProcessXML44::endXMLcommunityNode( ifstream& infile )

{
  string line;

  while( line.find( "</community>" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void ProcessXML44::endXMLtvegNode( ifstream& infile )

{
  string line;

  while( line.find( "</temveg>" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLcommunityNode( ifstream& infile,
                                       const string& rootnode )

{
  string line;

  int startString;

  string temp;

  string value;

  while( line.find( ">" ) == string::npos
          && !infile.eof() )
  {
    getline( infile, temp );

    if( temp.size() > 0 )
    {
      line += temp;
      temp.erase();
    }
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find community type in " << rootnode;
    cerr << " !!!" << endl;

    exit( -1 );
  }

  startString = line.find( "<community type = " );

  temp = line.substr( startString, 30 );

  startString = temp.find( '"' );

  value = temp.substr( (startString+1), 5 );
  
  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML44::getXMLcmntArrayDouble( ifstream& infile,
                                            const string& rootnode,
                                            const string& varnode,
                                            const int& index )

{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;


  while( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;

    exit( -1 );
  }

  startString = line.find( ">" );

  if( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );

    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLcmntArrayInt( ifstream& infile,
                                      const string& rootnode,
                                      const string& varnode,
                                      const int& index )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
          && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;

    exit( -1 );
  }

  startString = line.find( ">" );

  if( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );

    return atoi( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

long ProcessXML44::getXMLcmntArrayLong( ifstream& infile,
                                        const string& rootnode,
                                        const string& varnode,
                                        const int& index )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for community type ";
    cerr << index << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
  
    return atol( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML44::getXMLdouble( ifstream& infile,
                                   const string& rootnode,
                                   const string& varnode )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
    
    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLint( ifstream& infile,
                             const string& rootnode,
                             const string& varnode )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
  
    return atoi( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

long ProcessXML44::getXMLlong( ifstream& infile,
                               const string& rootnode,
                               const string& varnode )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode;
    cerr << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return -99; }
  else
  {
    value = line.substr( (startString+1), 20 );
  
    return atol( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLrootNode( ifstream& infile,
                                  const string& rootnode )

{
  string line;

  while( line.find( rootnode ) == string::npos && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << rootnode << " !!!" << endl;
    
    exit( -1 );
  }

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLsiteCommunityNode( ifstream& infile,
                                           const string& rootnode,
                                           string& description )

{
  int endString;

  string line;

  int startString;

  string temp;

  string temp2;

  string value;

  while( line.find( ">" ) == string::npos
         && !infile.eof() )
  {
    getline( infile, temp );
    
    if( temp.size() > 0 )
    {
      line += temp;
    
      temp.erase();
    }
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find community type in " << rootnode;
    cerr << " !!!" << endl;
    
    exit( -1 );
  }

  startString = line.find( "<community type = " );
  
  temp = line.substr( startString, 30 );
  
  startString = temp.find( '"' );
  
  value = temp.substr( (startString+1), 5 );

  temp.erase();
  
  startString = line.find( "description" );
  
  temp = line.substr( startString, 50 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 50 );
  
  endString = temp2.find( '"' );
  
  description = temp2.substr( 0, endString );

  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLsiteRootNode( ifstream& infile,
                                      const string& rootnode,
                                      string& version,
                                      string& sitename,
                                      string& sitecol,
                                      string& siterow,
                                      string& developer,
                                      string& updated )
{
  int endString;

  string line;

  int startString;

  string temp;

  string temp2;

  while( line.find( rootnode ) == string::npos && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << rootnode << " !!!" << endl;
    
    exit( -1 );
  }

  while( line.find( ">" ) == string::npos && !infile.eof() )
  {
    getline( infile, temp );
    
    if( temp.size() > 0 )
    {
      line += temp;
    
      temp.erase();
    }
  }


  startString = line.find( "version" );
  
  temp = line.substr( startString, 50 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 50 );
  
  endString = temp2.find( '"' );
  
  version = temp2.substr( 0, endString );

  
  temp.erase();
  
  startString = line.find( "site" );
  
  temp = line.substr( startString, 80 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 80 );
  
  endString = temp2.find( '"' );
  
  sitename = temp2.substr( 0, endString );

  
  temp.erase();
  
  startString = line.find( "longitude" );
  
  temp = line.substr( startString, 50 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 50 );
  
  endString = temp2.find( '"' );
  
  sitecol = temp2.substr( 0, endString );

  
  temp.erase();
  
  startString = line.find( "latitude" );
  
  temp = line.substr( startString, 50 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 50 );
  
  endString = temp2.find( '"' );
  
  siterow = temp2.substr( 0, endString );

  
  temp.erase();
  
  startString = line.find( "developedBy" );
  
  temp = line.substr( startString, 80 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 80 );
  
  endString = temp2.find( '"' );
  
  developer = temp2.substr( 0, endString );

  
  temp.erase();
  
  startString = line.find( "updated" );
  
  temp = line.substr( startString, 50 );
  
  startString = temp.find( '"' );
  
  temp2 = temp.substr( (startString+1), 50 );
  
  endString = temp2.find( '"' );
  
  updated = temp2.substr( 0, endString );

  return 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLtemvegNode( ifstream& infile,
                                    const string& rootnode )
{
  string line;
  
  int startString;

  string temp;

  string value;

  while( line.find( ">" ) == string::npos
         && !infile.eof() )
  {
    getline( infile, temp );
    
    if( temp.size() > 0 )
    {
      line += temp;
    
      temp.erase();
    }
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find TEMVEG type in " << rootnode;
    cerr << " !!!" << endl;
    
    exit( -1 );
  }

  startString = line.find( "<temveg type = " );
  
  temp = line.substr( startString, 30 );
  
  startString = temp.find( '"' );
  
  value = temp.substr( (startString+1), 5 );

  return atoi( value.c_str() );

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double ProcessXML44::getXMLtvegArrayDouble( ifstream& infile,
                                            const string& rootnode,
                                            const string& varnode,
                                            const int& index )

{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;


  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for TEMVEG = ";
    cerr << index << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
  
    return atof( value.c_str() );
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int ProcessXML44::getXMLtvegArrayInt( ifstream& infile,
                                      const string& rootnode,
                                      const string& varnode,
                                      const int& index )
{
  string endVarnode = "</" + varnode + ">";

  string line;

  int startString;

  string value;

  while( line.find( endVarnode ) == string::npos
         && !infile.eof() )
  {
    getline( infile, line );
  }

  if( 0 == line.size() )
  {
    cerr << endl << "Cannot find " << varnode << " for TEMVEG = ";
    cerr << index << " in " << rootnode << endl;
    
    exit( -1 );
  }

  startString = line.find( ">" );
  
  if( startString == string::npos ) { return (int) MISSING; }
  else
  {
    value = line.substr( (startString+1), 20 );
  
    return atoi( value.c_str() );
  }

};


