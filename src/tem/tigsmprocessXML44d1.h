/* *************************************************************
TIGSMPROCESSMXL44D1.H - Functions to process XML files for TEM

Programmer: Shaomin Hu and David Kicklighter
Creation Date: 17 June 2003

Modifications:

20060126 - DWK created by modifying ttprocessXML51.h
20060126 - DWK changed class ProcessXML50 to class ProcessXML43
20060126 - DWK changed include from temconsts51.hpp to 
           temconsts43.hpp
20080130 - DWK changed include from temconsts43.hpp to
           temigsmconsts44a.hpp
20080130 - DWK changed class ProcessXML43 to class ProcessXML44
20110707 - DWK changed include from temigsmconsts44a.hpp to
           temigsmconsts44c.hpp
20150429 - DWK changed include temigsmconsts44c.hpp to
            temigsmconstants.hpp
                        
************************************************************** */

#ifndef PROCESSXML44D_H
#define PROCESSXML44D_H

#include "temigsmconstants.hpp"

class ProcessXML44
{
  public:

  ProcessXML44();

  void endXMLcommunityNode( ifstream& infile );

  void endXMLtvegNode( ifstream& infile );

  double getXMLcmntArrayDouble( ifstream& infile,
                                const string& rootnode,
                                const string& varnode,
                                const int& index );

  int getXMLcmntArrayInt( ifstream& infile,
                          const string& rootnode,
                          const string& varnode,
                          const int& index );

  long getXMLcmntArrayLong( ifstream& infile,
                            const string& rootnode,
                            const string& varnode,
                            const int& index );

  int getXMLcommunityNode( ifstream& infile,
                           const string& rootnode );

  double getXMLdouble( ifstream& infile,
                       const string& rootnode,
                       const string& varnode );

  int getXMLint( ifstream& infile,
                 const string& rootnode,
                 const string& varnode );

  long getXMLlong( ifstream& infile,
                   const string& rootnode,
                   const string& varnode );

  int getXMLrootNode( ifstream& infile, const string& rootnode );

  int getXMLsiteCommunityNode( ifstream& infile,
                               const string& rootnode,
                               string& description );

  int getXMLsiteRootNode( ifstream& infile,
                          const string& rootnode,
                          string& version,
                          string& sitename,
                          string& sitecol,
                          string& siterow,
                          string& developer,
                          string& updated );

  int getXMLtemvegNode( ifstream& infile, const string& rootnode );

  double getXMLtvegArrayDouble( ifstream& infile,
                                const string& rootnode,
                                const string& varnode,
                                const int& index );

  int getXMLtvegArrayInt( ifstream& infile,
                          const string& rootnode,
                          const string& varnode,
                          const int& index );

};

#endif
