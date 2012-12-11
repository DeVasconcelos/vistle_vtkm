#ifndef READVISTLE_H
#define READVISTLE_H

#include <string.h>
#include "module.h"

class ReadVistle: public vistle::Module {

 public:
   ReadVistle(const std::string &shmname, int rank, int size, int moduleID);
   ~ReadVistle();

 private:
   bool load(const std::string & name);
   virtual bool compute();
   bool first_object;
};

#endif
