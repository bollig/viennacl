#include <viennacl/distributed/devices_grid.hpp>


int main(){
    viennacl::distributed::devices_grid grid = viennacl::distributed::current_default_grid();
    std::cout << grid.size1() << " " << grid.size2() << std::endl;
}
