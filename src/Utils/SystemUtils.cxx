//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 08, 2009 - CA
   That file was added in 2.5.1

*/
//____________________________________________________________________________

#include <sys/types.h>
#include <dirent.h>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"

//___________________________________________________________________________
vector<string> 
  genie::utils::system::GetAllFilesInPath(string path, string extension)
{
  vector<string> files;

  DIR *dp = NULL;
  struct dirent *dirp = NULL;
  if((dp = opendir(path.c_str())) == NULL) {
    LOG("System", pERROR) << "Can not open directory: " << path;
  }
  else {
    while ((dirp = readdir(dp)) != NULL) {
      string filename = path + "/" + string(dirp->d_name);
      bool match = false;
      if(extension.size()==0) match=true; // get all files
      else {
        // match extension
        vector<string> filenamev = genie::utils::str::Split(filename,".");
        if(filenamev.size()==0) break;
        string file_extension = filenamev[filenamev.size()-1];
        match = (file_extension.find(extension) != string::npos);
      }
      if(match) {
        files.push_back(filename);
      }
    }
    closedir(dp);
  }
  return files;
}
//___________________________________________________________________________