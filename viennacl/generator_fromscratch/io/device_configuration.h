#ifndef VIENNACL_GENERATOR_DEVICE_CONFIGURATION_H
#define VIENNACL_GENERATOR_DEVICE_CONFIGURATION_H

#include "pugixml/src/pugixml.hpp"


namespace viennacl{

namespace generator{

namespace io
{
  namespace tag
  {
    static std::string root     = "device_configuration";
    static std::string devices  = "devices";
    static std::string device   = "device";
    static std::string name     = "name";
    static std::string driver   = "driver";
    static std::string kernel_type    = "kernel_type";
  }

  namespace val
  {
      static std::string inner_product = "inner_product";
      static std::string globsize = "globalsize";
      static std::string locsize  = "localsize";

  }

  } // end namespace io


}

}
#endif // DEVICE_CONFIGURATION_H
