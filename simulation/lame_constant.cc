#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <iostream>
#include <map>
#include <string>

using namespace dealii;

struct MaterialProperty
{
  double rho; // [g/cm^3]
  double vp;  // [m/s]
  double vs;  // [m/s]
};

template <int dim>
class ElasticWave
{
public:
  ElasticWave() = default;
  void run();

private:
  void compute_lame_constants(const MaterialProperty &mat,
                             double &lambda,
                             double &mu);

  void print_material_lame();
};

template <int dim>
void ElasticWave<dim>::compute_lame_constants(const MaterialProperty &mat,
                                             double &lambda,
                                             double &mu)
{
  const double rho_SI = mat.rho * 1000.0; // [kg/m^3]
  mu = rho_SI * mat.vs * mat.vs;
  lambda = rho_SI * mat.vp * mat.vp - 2.0 * mu;
}

template <int dim>
void ElasticWave<dim>::print_material_lame()
{
  std::map<std::string, MaterialProperty> materials;

  // 氷層
  materials["ice"] = {0.9, 3600.0, 1800.0}; // 中間値
  // 両脇層
  materials["saturated_shale_and_silt"] = {2.25, 1850.0, 625.0}; // 中間値
  // 下部層
  materials["gneiss"] = {2.6, 5250.0, 2950.0}; // 中間値

  for (const auto &entry : materials)
  {
    const std::string &name = entry.first;
    const MaterialProperty &mat = entry.second;
    double lambda = 0.0, mu = 0.0;

    compute_lame_constants(mat, lambda, mu);

    std::cout << "Material: " << name << std::endl;
    std::cout << "  rho [g/cm^3]: " << mat.rho << std::endl;
    std::cout << "  Vp  [m/s]:    " << mat.vp << std::endl;
    std::cout << "  Vs  [m/s]:    " << mat.vs << std::endl;
    std::cout << "  lambda [Pa]:  " << lambda << std::endl;
    std::cout << "  mu     [Pa]:  " << mu << std::endl;
    std::cout << std::endl;
  }
}

template <int dim>
void ElasticWave<dim>::run()
{
  deallog.depth_console(0);
  print_material_lame();
}

int main()
{
  try
  {
    constexpr int dim = 2;
    ElasticWave<dim> elastic_problem;
    elastic_problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << "Exception: " << exc.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception!" << std::endl;
    return 1;
  }

  return 0;
}
