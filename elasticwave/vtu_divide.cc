//上下でλ、μの値を変えてみた。
//上のほうが値が大きくなるように調整して上側の波の伝搬のほうが早くなるようにしてある

//有限要素なのに波源を１点に指定したために不安定になっており波が安定していない

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// λとμを空間依存にする関数クラス
template <int dim>
class LambdaMuFunction : public Function<dim>
{
public:
  LambdaMuFunction() : Function<dim>(2) {} // 2成分: 0番目がlambda, 1番目がmu

  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2));
    if (p[1] < 0.5)
    {
      values[0] = 6.0; // lambda 下部
      values[1] = 3.0; // mu 下部
    }
    else
    {
      values[0] = 6.0; // lambda 上部
      values[1] = 3.0; // mu 上部
    }
  }
};

template <int dim>
class ElasticWaveProblem
{
public:
  ElasticWaveProblem(const unsigned int degree);
  void run();

private:
  void setup_system();
  void assemble_system();
  void time_step();
  void output_results(const unsigned int timestep) const;

  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern  sparsity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> stiffness_matrix;

  Vector<double> displacement;
  Vector<double> velocity;
  Vector<double> acceleration;
  Vector<double> system_rhs;

  double time;
  double time_step_size;
  unsigned int timestep_number;
  const unsigned int n_timesteps;

  const double rho = 1.0; // 密度は均一と仮定

  LambdaMuFunction<dim> lambda_mu; // 空間依存λ,μを管理する関数
};

template <int dim>
ElasticWaveProblem<dim>::ElasticWaveProblem(const unsigned int degree)
  :
  fe(FE_Q<dim>(degree), dim),
  dof_handler(triangulation),
  time(0.0),
  time_step_size(1e-3),
  timestep_number(0),
  n_timesteps(800)
{}

template <int dim>
void ElasticWaveProblem<dim>::setup_system()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(6);

  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);
  stiffness_matrix.reinit(sparsity_pattern);

  displacement.reinit(dof_handler.n_dofs());
  velocity.reinit(dof_handler.n_dofs());
  acceleration.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void ElasticWaveProblem<dim>::assemble_system()
{
  mass_matrix = 0;
  stiffness_matrix = 0;

  QGauss<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_JxW_values | update_quadrature_points);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  Vector<double> lambda_mu_values(2);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    cell_mass_matrix = 0;
    cell_stiffness_matrix = 0;

    fe_values.reinit(cell);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const Point<dim> &x_q = fe_values.quadrature_point(q);

      lambda_mu.vector_value(x_q, lambda_mu_values);
      // const double lambda = lambda_mu_values[0];
      const double mu = lambda_mu_values[1];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          const unsigned int component_j = fe.system_to_component_index(j).first;

          // 質量行列 (ρ * φ_i * φ_j)
          if (component_i == component_j)
          {
            cell_mass_matrix(i, j) +=
                (rho * fe_values.shape_value(i, q) *
                 fe_values.shape_value(j, q) * fe_values.JxW(q));
          }

          // 剛性行列（簡易版）
          cell_stiffness_matrix(i, j) +=
              (mu * fe_values.shape_grad(i, q)[component_i] *
               fe_values.shape_grad(j, q)[component_j] * fe_values.JxW(q));
        }
      }
    }

    cell->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        mass_matrix.add(local_dof_indices[i], local_dof_indices[j],
                        cell_mass_matrix(i, j));
        stiffness_matrix.add(local_dof_indices[i], local_dof_indices[j],
                             cell_stiffness_matrix(i, j));
      }
    }
  }
}

template <int dim>
void ElasticWaveProblem<dim>::time_step()
{
  Vector<double> tmp(stiffness_matrix.m());
  stiffness_matrix.vmult(tmp, displacement);

  SparseDirectUMFPACK solver;
  solver.initialize(mass_matrix);
  solver.vmult(acceleration, tmp);
  acceleration *= -1;

  Vector<double> displacement_new(displacement.size());
  Vector<double> velocity_new(velocity.size());

  for (unsigned int i = 0; i < displacement.size(); ++i)
  {
    displacement_new[i] = displacement[i] + time_step_size * velocity[i] +
                          0.5 * time_step_size * time_step_size * acceleration[i];
    velocity_new[i] = velocity[i] + 0.5 * time_step_size * acceleration[i];
  }

  stiffness_matrix.vmult(tmp, displacement_new);
  solver.vmult(acceleration, tmp);
  acceleration *= -1;

  for (unsigned int i = 0; i < velocity.size(); ++i)
  {
    velocity_new[i] += 0.5 * time_step_size * acceleration[i];
  }

  displacement = displacement_new;
  velocity = velocity_new;

  time += time_step_size;
  ++timestep_number;
}

template <int dim>
void ElasticWaveProblem<dim>::output_results(const unsigned int timestep) const
{
  DataOut<dim> data_out;

  std::vector<std::string> solution_names(dim, "displacement");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(displacement, solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);

  data_out.build_patches();

  const std::string filename =
      "solution-" + Utilities::int_to_string(timestep, 4) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
  // ファイル名をコンソールに出力
  std::cout << "Output written to file: " << filename << std::endl;
}

template <int dim>
void ElasticWaveProblem<dim>::run()
{
  setup_system();
  assemble_system();

  // 初期条件
  displacement = 0;
  velocity = 0;
  acceleration = 0;

  // 波源を y=0.5 の境界付近に設定
  const Point<dim> source_point(0, 0.5);

  std::vector<Point<dim>> support_points(dof_handler.n_dofs());
  MappingQ1<dim> mapping;
  DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

  double min_dist = 1e10;
  unsigned int closest_dof = 0;

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
  {
    double dist = support_points[i].distance(source_point);
    if (dist < min_dist)
    {
      min_dist = dist;
      closest_dof = i;
    }
  }
/////////////////////////////////////////////////////////////////////////
  // displacement[closest_dof] = 1.0e-3; // 小さな初期変位（波源）

  // for (unsigned int step = 0; step < n_timesteps; ++step)
  // {
  //   time_step();

  //   if (step % 10 == 0)
  //     output_results(step);
  // }
// ///////////////////////////////////////////////////////////////////////////
  const double source_amplitude = 1.0e-3;
  const double frequency = 50.0; // Hz
  const double omega = 2.0 * numbers::PI * frequency;

  double T_pulse = 5.0 / frequency; // 周波数からパルス長を決める

  for (unsigned int step = 0; step < n_timesteps; ++step)
  {
      if (time <= T_pulse) // 一定時間だけ強制
      {
          displacement[closest_dof] = source_amplitude * std::sin(omega * time);
          velocity[closest_dof] = omega * source_amplitude * std::cos(omega * time);
      }
      else
      {
          // 波源拘束を外す
          displacement[closest_dof] = 0.0;
          velocity[closest_dof] = 0.0;
      }

      time_step();

      if (step % 10 == 0)
          output_results(step);
  }

/////////////////////////////////////////////////////////////////////////////////


}

int main()
{
  try
  {
    deallog.depth_console(2);

    ElasticWaveProblem<2> elastic_wave_problem(1);
    elastic_wave_problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << "\n\nException on processing: " << exc.what() << "\n";
    return 1;
  }
  catch (...)
  {
    std::cerr << "\n\nUnknown exception!\n";
    return 1;
  }

  return 0;
}
