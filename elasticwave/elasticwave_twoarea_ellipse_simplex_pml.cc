#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

namespace ElasticWave
{
using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide() : Function<dim>(dim) {}

  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    values = 0;
    const double t = this->get_time();

    const double pulse_amplitude = 10.0;
    const double pulse_duration = 0.1;
    const double r = p.norm();

    if (r < 0.1)
    {
      const double arg = numbers::PI * (t - pulse_duration) / pulse_duration;
      const double f = pulse_amplitude * std::exp(-arg * arg);
      values[1] = f;
    }
  }
};

template <int dim>
struct MaterialProperty
{
  static double rho(const types::material_id id)
  {
    if (id == 0) return 1.0;  // 波源
    if (id == 1) return 10.0; // 構造体
    if (id == 2) return 1.0;  // PML
    return 1.0;
  }

  static double lambda(const types::material_id id)
  {
    if (id == 0) return 1.0;
    if (id == 1) return 10.0;
    if (id == 2) return 1.0;
    return 1.0;
  }

  static double mu(const types::material_id id)
  {
    if (id == 0) return 1.0;
    if (id == 1) return 10.0;
    if (id == 2) return 1.0;
    return 1.0;
  }
};

template <int dim>
class ElasticWaveProblem
{
public:
  ElasticWaveProblem();
  void run();

private:
  void setup_system();
  void assemble_system();
  void solve_timestep();
  void output_results(const unsigned int timestep_number) const;

  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> laplace_matrix;

  Vector<double> solution;
  Vector<double> old_solution;
  Vector<double> old_old_solution;
  Vector<double> system_rhs;

  double time;
  double time_step;
  unsigned int timestep_number;
  const unsigned int n_cycles = 200;
};

template <int dim>
ElasticWaveProblem<dim>::ElasticWaveProblem()
  : fe(FE_Q<dim>(1), dim), dof_handler(triangulation), time(0), time_step(1e-2), timestep_number(0)
{}

template <int dim>
void ElasticWaveProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);
  laplace_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  old_old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void ElasticWaveProblem<dim>::assemble_system()
{
  mass_matrix = 0;
  laplace_matrix = 0;

  QGauss<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients |
                          update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_laplace_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_mass_matrix = 0;
    cell_laplace_matrix = 0;

    const auto mat_id = cell->material_id();
    const double rho = MaterialProperty<dim>::rho(mat_id);
    const double lambda = MaterialProperty<dim>::lambda(mat_id);
    const double mu = MaterialProperty<dim>::mu(mat_id);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int comp_i = fe.system_to_component_index(i).first;
        const Tensor<1, dim> phi_i = fe_values[FEValuesExtractors::Vector(0)].value(i, q);
        const Tensor<2, dim> grad_phi_i = fe_values[FEValuesExtractors::Vector(0)].gradient(i, q);

        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          const unsigned int comp_j = fe.system_to_component_index(j).first;
          const Tensor<1, dim> phi_j = fe_values[FEValuesExtractors::Vector(0)].value(j, q);
          const Tensor<2, dim> grad_phi_j = fe_values[FEValuesExtractors::Vector(0)].gradient(j, q);

          double mass = rho * phi_i[comp_i] * phi_j[comp_j] * fe_values.JxW(q);
          cell_mass_matrix(i, j) += mass;

          double laplace = (lambda * trace(grad_phi_i) * trace(grad_phi_j) +
                            2.0 * mu * scalar_product(grad_phi_i, grad_phi_j)) * fe_values.JxW(q);
          cell_laplace_matrix(i, j) += laplace;
        }
      }
    }

    cell->get_dof_indices(local_dof_indices);
    mass_matrix.add(local_dof_indices, cell_mass_matrix);
    laplace_matrix.add(local_dof_indices, cell_laplace_matrix);
  }
}

template <int dim>
void ElasticWaveProblem<dim>::solve_timestep()
{
  RightHandSide<dim> rhs_function;
  rhs_function.set_time(time);

  system_rhs = 0;

  Vector<double> tmp(dof_handler.n_dofs());
  VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                      rhs_function, tmp);
  system_rhs = tmp;

  Vector<double> acceleration(dof_handler.n_dofs());

  mass_matrix.vmult(acceleration, solution);
  laplace_matrix.vmult_add(acceleration, solution);
  acceleration.add(1.0, system_rhs);

  Vector<double> new_solution(dof_handler.n_dofs());
  new_solution = 0;

  for (unsigned int i = 0; i < new_solution.size(); ++i)
  {
    new_solution[i] = (2.0 * old_solution[i] - old_old_solution[i] +
                       time_step * time_step * acceleration[i]);
  }

  old_old_solution = old_solution;
  old_solution = new_solution;
  solution = new_solution;
}

template <int dim>
void ElasticWaveProblem<dim>::output_results(const unsigned int timestep_number) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  std::vector<std::string> solution_names(dim, "displacement");
  data_out.add_data_vector(solution, solution_names);
  data_out.build_patches();
  std::ofstream output("solution-" + std::to_string(timestep_number) + ".vtu");
  data_out.write_vtu(output);
}

template <int dim>
void ElasticWaveProblem<dim>::run()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);

  for (const auto &cell : triangulation.active_cell_iterators())
  {
    const Point<dim> p = cell->center();
    const double x = p[0];
    const double y = p[1];
    const double a = 0.5;
    const double b = 0.3;

    if ((x * x) / (a * a) + (y * y) / (b * b) < 1.0)
      cell->set_material_id(0); // 中央波源
    else if (std::abs(x) >= 0.9 || std::abs(y) >= 0.9)
      cell->set_material_id(2); // PML
    else
      cell->set_material_id(1); // 構造体
  }

  setup_system();
  assemble_system();

  output_results(timestep_number);

  for (timestep_number = 1; timestep_number <= n_cycles; ++timestep_number)
  {
    time += time_step;
    solve_timestep();
    output_results(timestep_number);
  }
}

} // namespace ElasticWave

int main()
{
  try
  {
    using namespace ElasticWave;
    ElasticWaveProblem<2> wave_problem;
    wave_problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << "Exception: " << exc.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << "Unknown exception!" << std::endl;
    return 1;
  }

  return 0;
}
