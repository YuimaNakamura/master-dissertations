// ヘッダー部（元コードと同じ）
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

namespace ElasticWave2D
{
  using namespace dealii;

  template <int dim>
  struct MaterialProperty
  {
    static double rho(const types::material_id id)
    {
      return (id == 0) ? 1.0 : 1000.0;
    }

    static double lambda(const types::material_id id)
    {
      return (id == 0) ? 1.0 : 10.0;
    }

    static double mu(const types::material_id id)
    {
      return (id == 0) ? 1.0 : 10.0;
    }
  };

  template <int dim>
  class ElasticWave
  {
  public:
    ElasticWave();
    void run();

  private:
    void setup_system();
    void assemble_system();
    void assemble_mass_matrix();
    void initialize_solution();
    void time_step();
    void output_results(const unsigned int timestep) const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FESystem<dim>      fe;
    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    SparseMatrix<double> mass_matrix;

    Vector<double> solution_n, solution_nm1, solution_np1;
    Vector<double> system_rhs;
    Vector<double> mass_matrix_diagonal_inverse;

    double time;
    double time_step_size;
    unsigned int timestep_number;
  };

  template <int dim>
  ElasticWave<dim>::ElasticWave()
    : dof_handler(triangulation)
    , fe(FE_Q<dim>(1), dim)
    , time(0.)
    , time_step_size(1e-3)
    , timestep_number(0)
  {}

  template <int dim>
  void ElasticWave<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution_n.reinit(dof_handler.n_dofs());
    solution_nm1.reinit(dof_handler.n_dofs());
    solution_np1.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(dim), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    stiffness_matrix.reinit(sparsity_pattern);
    mass_matrix.reinit(sparsity_pattern);
  }

  template <int dim>
  void ElasticWave<dim>::assemble_system()
  {
    stiffness_matrix = 0;
    QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points = quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      const types::material_id id = cell->material_id();

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double lambda = MaterialProperty<dim>::lambda(id);
        const double mu     = MaterialProperty<dim>::mu(id);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int ci = fe.system_to_component_index(i).first;
          const Tensor<1, dim> grad_i = fe_values.shape_grad(i, q);

          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            const unsigned int cj = fe.system_to_component_index(j).first;
            const Tensor<1, dim> grad_j = fe_values.shape_grad(j, q);

            double val = lambda * grad_i[ci] * grad_j[cj]
                       + mu * (grad_i * grad_j) * (ci == cj)
                       + mu * grad_i[cj] * grad_j[ci];

            cell_matrix(i, j) += val * fe_values.JxW(q);
          }
        }
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, local_dof_indices, stiffness_matrix);
    }
  }

  template <int dim>
  void ElasticWave<dim>::assemble_mass_matrix()
  {
    mass_matrix = 0;
    QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points = quadrature_formula.size();
    FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_mass = 0;
      const types::material_id id = cell->material_id();

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double rho = MaterialProperty<dim>::rho(id);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            cell_mass(i, j) += rho * fe_values.shape_value(i, q) *
                               fe_values.shape_value(j, q) *
                               fe_values.JxW(q);
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_mass, local_dof_indices, mass_matrix);
    }

    mass_matrix_diagonal_inverse.reinit(mass_matrix.m());
    for (unsigned int i = 0; i < mass_matrix.m(); ++i)
      mass_matrix_diagonal_inverse[i] = 1.0 / mass_matrix.diag_element(i);
  }

  template <int dim>
  void ElasticWave<dim>::initialize_solution()
  {
    solution_n = 0;
    solution_nm1 = 0;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        const unsigned int comp = fe.system_to_component_index(i).first;
        const Point<dim> p = cell->vertex(i / dim); // 簡易な近似

        if (p.distance(Point<dim>()) < 0.2 && comp == 0)
          solution_n(local_dof_indices[i]) = 0.01 * std::exp(-50 * p.square());
      }
    }

    solution_nm1 = solution_n;
  }

  template <int dim>
  void ElasticWave<dim>::time_step()
  {
    stiffness_matrix.vmult(system_rhs, solution_n);
    system_rhs *= -1.0;

    for (unsigned int i = 0; i < system_rhs.size(); ++i)
      solution_np1[i] = mass_matrix_diagonal_inverse[i] * system_rhs[i] * time_step_size * time_step_size
                        + 2.0 * solution_n[i] - solution_nm1[i];

    constraints.distribute(solution_np1);
    solution_nm1 = solution_n;
    solution_n = solution_np1;
    time += time_step_size;
    ++timestep_number;
  }

  template <int dim>
  void ElasticWave<dim>::output_results(const unsigned int timestep) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_n, "displacement");
    data_out.build_patches();

    const std::string filename = "solution-" + std::to_string(timestep) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);

    std::cout << "VTU file written: " << filename << std::endl;
  }

  template <int dim>
  void ElasticWave<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(6);

    // material_id の設定（波打つ境界 y < 0.3 sin(3πx)）
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      const Point<dim> p = cell->center();
      if (p[1] < 0.3 * std::sin(3 * numbers::PI * p[0]))
        cell->set_material_id(0); // 柔らかい
      else
        cell->set_material_id(1); // 硬い
    }

    setup_system();
    assemble_system();
    assemble_mass_matrix();
    initialize_solution();
    output_results(timestep_number);

    const unsigned int n_time_steps = 1000;
    for (unsigned int step = 1; step <= n_time_steps; ++step)
    {
      time_step();
      if (step % 10 == 0)
        output_results(timestep_number);
    }
  }
}

int main()
{
  try
  {
    ElasticWave2D::ElasticWave<2> simulation;
    simulation.run();
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
