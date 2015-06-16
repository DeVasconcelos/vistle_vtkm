#include <sstream>
#include <iomanip>

#include <core/object.h>
#include <core/geometry.h>
#include <core/normals.h>
#include <core/texture1d.h>

#include "Collect.h"

MODULE_MAIN(Collect)

using namespace vistle;


Collect::Collect(const std::string &shmname, const std::string &name, int moduleID)
   : Module("Collect", shmname, name, moduleID) {

   createInputPort("grid_in");
   createInputPort("normal_in");
   createInputPort("texture_in");

   createOutputPort("grid_out");
}

Collect::~Collect() {

}


bool Collect::compute() {

   vistle::Object::const_ptr grid = expect<Object>("grid_in");
   if (!grid)
      return true;

   vistle::Geometry::ptr geom(new vistle::Geometry(grid));
   geom->setMeta(grid->meta());

   vistle::Object::const_ptr norm = accept<Normals>("normal_in");
   if (norm)
      geom->setNormals(norm);

   vistle::Object::const_ptr tex = accept<Texture1D>("texture_in");
   if (tex)
      geom->setTexture(tex);

   addObject("grid_out", geom);

   return true;
}
