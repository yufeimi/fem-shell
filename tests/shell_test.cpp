#include "fem-shell.h"

// Begin the main program.
int main(int argc, char **argv)
{
  shellparam param;
  // read command-line arguments and initialize global variables
  if (param.read_parameters(argc, argv))
    {
      std::cout << "Read command-line arguments.......OK" << std::endl;
    }
  else
    {
      std::cout << "Read command-line arguments.......FAILED" << std::endl;
      return -1;
    }

  // Initialize libMesh and any dependent library
  LibMeshInit init(argc, argv);

  // Skip this program if libMesh was compiled as 1D-only.
  libmesh_example_requires(LIBMESH_DIM >= 2, "2D support");

  // Initialize the mesh
  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), 2);
  mesh.allow_renumbering(
    false); // prevent libMesh from renumber nodes on its own
  mesh.read(param.in_filename);
  // throws error if file does not exist and checks content
  // Print information about the mesh to the screen
  mesh.print_info();

  shellsolid shell(mesh, param);
  shell.run();
  return 0;
}
