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

  Triangulation<dim>         triangulation;
  FESystem<dim>              fe;
  DoFHandler<dim>            dof_handler;

  SparsityPattern            sparsity_pattern;
  SparseMatrix<double>       mass_matrix;
  SparseMatrix<double>       stiffness_matrix;

  Vector<double>             displacement;
  Vector<double>             velocity;
  Vector<double>             acceleration;
  Vector<double>             system_rhs;

  double                     time;
  double                     time_step_size;
  unsigned int               timestep_number;
  const unsigned int         n_timesteps;

  // 物理パラメータ
  double                     rho;    // 密度
  double                     lambda; // ラメ定数λ
  double                     mu;     // ラメ定数μ
};


template <int dim>
ElasticWaveProblem<dim>::ElasticWaveProblem(const unsigned int degree)
  :
  fe(FE_Q<dim>(degree), dim),
  dof_handler(triangulation),
  time(0.0),
  time_step_size(1e-3),
  timestep_number(0),
  n_timesteps(500),
  rho(1.0),
  lambda(1.0),
  mu(1.0)
{}


template <int dim>
void ElasticWaveProblem<dim>::setup_system()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(5);

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
                          update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    cell_mass_matrix = 0;
    cell_stiffness_matrix = 0;

    fe_values.reinit(cell);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;

        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          const unsigned int component_j = fe.system_to_component_index(j).first;

          // 質量行列 (ρ * φ_i * φ_j) を対角成分で組み立て
          if (component_i == component_j)
            cell_mass_matrix(i, j) +=
                (rho * fe_values.shape_value(i, q) *
                 fe_values.shape_value(j, q) * fe_values.JxW(q));

          // 剛性行列の簡易版 (λ, μを使った応力-ひずみ関係はここでは省略して簡単化)
          // ここでは単純なラプラシアン近似での剛性行列に置き換え（本格実装は step-20, step-62 参照）
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
}

template <int dim>
void ElasticWaveProblem<dim>::run()
{
  setup_system();
  assemble_system();

  // 初期条件：全て0ではなく、中央付近の自由度に小さな変位を設定する
  displacement = 0;
  velocity = 0;
  acceleration = 0;

  // 中央点の座標 (0.5, 0.5) に最も近い自由度にだけ値を入れる例（2Dの場合）
  const Point<dim> center(0.5, 0.5);

  for (unsigned int i = 0; i < dof_handler.n_dofs(); i += dim)
  {
    // 各自由度の空間座標を取得
    const auto support_points = dof_handler.get_fe().get_unit_support_points();
    // support_pointsは単位セルなのでグローバル座標取得用にCellでやる必要あり
    // 簡単化のためここでは近似的にDoFの位置は計算しない。正確にしたい場合はCellを経由してください

    // 代わりにDoFの座標はdof_handlerで取得可能（以下参照）
  }

  // 以下は正確にDoF座標を取得して中央に近い自由度を探すサンプル
  std::vector<Point<dim>> support_points(dof_handler.n_dofs());
  MappingQ1<dim> mapping;
  DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

  // 中央に最も近い自由度を探し、そこに変位値を与える（x成分のみ例）
  double min_dist = 1e10;
  unsigned int closest_dof = 0;

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
  {
    const double dist = support_points[i].distance(center);
    if (dist < min_dist)
    {
      min_dist = dist;
      closest_dof = i;
    }
  }

  displacement[closest_dof] = 1.0e-3; // 0.001だけ変位を与える（例）

  for (unsigned int step = 0; step < n_timesteps; ++step)
  {
    time_step();

    if (step % 10 == 0)
      output_results(step);
  }
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
