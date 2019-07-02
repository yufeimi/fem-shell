#include "fem-shell.h"

// Begin the main program.
int main(int argc, char **argv)
{
  using namespace ShellSolid;
  using namespace libMesh;

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
  // Hard coded parameters for test C
  param.nu = 0.3;
  param.em = 10.92;
  param.thickness = 1.0;
  param.debug = false;
  std::cout << "Outfile: " << param.out_filename << " " << param.isOutfileSet
            << std::endl;
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

  // Construct BC map for displacement
  std::map<boundary_id_type, unsigned int> dirichlet_bcs;
  // 1-x, 2-y, 3-xy, 4-z, 5-xz, 6-yz, 7-xyz
  dirichlet_bcs.insert(std::pair<boundary_id_type, unsigned int>(0, 4));
  shellsolid shell(mesh, param);
  shell.make_constraints(dirichlet_bcs);
  shell.run();
  return 0;
}
