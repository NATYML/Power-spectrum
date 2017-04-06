#include "/home/nataly/Programs/FOFReaderLib-master/FOFReaderLib/FOFReaderLib.h"

int main(void)
{
    // Simple line to read a Field. The second arg is to ask reading the ids (default is false)
    //FOFCube cube("/path_for_data/fof_boxlen648_n1024_lcdmw5_cube_00001", true);
    FOFCube cube("./fof_boxlen162_n1024_lcdmw5_cube_00000", true);
    

    // Show some details about the read cube as an example of usage
    std::cout
    << cube.npart() << " particles, "
    << "area: (" << cube.minX() << "," << cube.minY() << "," << cube.minZ() << ") "
    << "to (" << cube.maxX() << "," << cube.maxY() << "," << cube.maxZ() << ")"
    << std::endl;
   for(int j=0; j< std::max(10,cube.npart()); j++) {
    std::cout <<   j << " "
    //<< "id: " << cube.id(j) << " "
    <<  cube.posX(j) << " " << cube.posY(j) << " " << cube.posZ(j) 
    //<< "velocity (" << cube.velX(j) << "," << cube.velY(j) << "," << cube.velZ(j) << ")"
    << std::endl;
    }
    
return 0;
}