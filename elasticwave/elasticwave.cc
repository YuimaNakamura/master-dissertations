/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2000 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>

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

    FESystem<dim> fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    SparseMatrix<double> mass_matrix;

    Vector<double> solution_n;    // u^n (変位 現在時刻)
    Vector<double> solution_nm1;  // u^{n-1} (一つ前の時刻)
    Vector<double> solution_np1;  // u^{n+1} (次の時刻)

    Vector<double> system_rhs;

    Vector<double> mass_matrix_diagonal_inverse;

    double time;
    double time_step_size;
    unsigned int timestep_number;

    // 物理パラメータ（任意に調整可能）
    const double rho;      // 密度
    const double lambda;   // ラメ定数 λ
    const double mu;       // ラメ定数 μ
  };


  template <int dim>
  ElasticWave<dim>::ElasticWave()
    : dof_handler(triangulation)
    , fe(FE_Q<dim>(1), dim) // 2次元なら2成分分の1次要素
    , time(0.)
    , time_step_size(1e-3) // 時間刻み（要調整）
    , timestep_number(0)
    , rho(1.0)
    , lambda(1.0)
    , mu(1.0)
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

    // 境界ID=0の境界でゼロ変位（固定）
    VectorTools::interpolate_boundary_values(
      dof_handler,
      0,
      Functions::ZeroFunction<dim>(dim),
      constraints);

    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    stiffness_matrix.reinit(sparsity_pattern);
    mass_matrix.reinit(sparsity_pattern);
  }

  // 剛性行列の組み立て（線形弾性）
  template <int dim>
  void ElasticWave<dim>::assemble_system()
  {
    stiffness_matrix = 0;

    QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const unsigned int component_i =
                  fe.system_to_component_index(i).first;

                const Tensor<1, dim> phi_i_grad =
                  fe_values.shape_grad(i, q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    const unsigned int component_j =
                      fe.system_to_component_index(j).first;

                    const Tensor<1, dim> phi_j_grad =
                      fe_values.shape_grad(j, q);

                    double val = 0.0;

                    // 線形弾性の剛性行列の積分
                    val +=
                      (lambda * phi_i_grad[component_i] * phi_j_grad[component_j]
                       + mu * (phi_i_grad * phi_j_grad) * ((component_i == component_j) ? 1.0 : 0.0)
                       + mu * phi_i_grad[component_j] * phi_j_grad[component_i]);

                    cell_matrix(i, j) += val * fe_values.JxW(q);
                  }
              }
          }
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, local_dof_indices, stiffness_matrix);
      }
  }

  // 質量行列組み立て（対角行列近似）
  template <int dim>
  void ElasticWave<dim>::assemble_mass_matrix()
  {
    mass_matrix = 0;

    QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_mass = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    // M_ij = \int rho * N_i * N_j dx
                    cell_mass(i, j) +=
                      rho * fe_values.shape_value(i, q) *
                      fe_values.shape_value(j, q) *
                      fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_mass, local_dof_indices, mass_matrix);
      }

    // 質量行列の対角成分の逆数を計算（陽的時間積分用）
    mass_matrix_diagonal_inverse.reinit(mass_matrix.m());

    for (unsigned int i = 0; i < mass_matrix.m(); ++i)
      {
        const double diag = mass_matrix.diag_element(i);
        Assert(diag > 1e-15, ExcMessage("質量行列の対角成分が小さすぎます"));
        mass_matrix_diagonal_inverse[i] = 1.0 / diag;
      }
  }

  // 初期条件：ここではゼロ変位＋簡単な局所パルスを初期変位に設定する例
  template <int dim>
  void ElasticWave<dim>::initialize_solution()
  {
    solution_n = 0;
    solution_nm1 = 0;
    solution_np1 = 0;

    // 例: 中央付近に小さい初期変位パルスを与える（x座標で）
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        std::vector<Point<dim>> support_points(fe.n_dofs_per_cell());
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          support_points[i] = cell->vertex(i / dim); // 簡易的に頂点座標を取得

        std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            const Point<dim> &p = support_points[i];

            // 中央付近 (x,y)=(0,0) に局所パルス
            if (p.distance(Point<dim>()) < 0.2)
              {
                if (component_i == 0)
                  solution_n(local_dof_indices[i]) = 0.01 * std::exp(-50 * p.square());
                else
                  solution_n(local_dof_indices[i]) = 0.0;
              }
          }
      }

    // 1つ前の時刻は0で初期化（静止状態）
    solution_nm1 = solution_n;
  }

  template <int dim>
  void ElasticWave<dim>::time_step()
  {
    // system_rhs = - K * u^n (外力は0とする)
    stiffness_matrix.vmult(system_rhs, solution_n);
    system_rhs *= -1.0;

    // 陽的中心差分法による時間積分
    for (unsigned int i = 0; i < system_rhs.size(); ++i)
      {
        solution_np1[i] = mass_matrix_diagonal_inverse[i] * system_rhs[i] * time_step_size * time_step_size
                          + 2.0 * solution_n[i]
                          - solution_nm1[i];
      }

    // 境界条件を満たすように強制的にセット（ゼロ固定境界）
    constraints.distribute(solution_np1);

    // 時刻更新
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

    // 変位ベクトルはdim成分のベクトルなので一括で追加
    data_out.add_data_vector(solution_n, "displacement");

    data_out.build_patches();

    std::string filename = "solution-" + std::to_string(timestep) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
    }


  template <int dim>
  void ElasticWave<dim>::run()
  {
    std::cout << "Generating mesh..." << std::endl;
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(6); // メッシュ細かさは要調整

    std::cout << "Setting up system..." << std::endl;
    setup_system();

    std::cout << "Assembling stiffness matrix..." << std::endl;
    assemble_system();

    std::cout << "Assembling mass matrix..." << std::endl;
    assemble_mass_matrix();

    std::cout << "Initializing solution..." << std::endl;
    initialize_solution();

    output_results(timestep_number);

    const unsigned int n_time_steps = 500; // 総時間ステップ数

    std::cout << "Starting time stepping..." << std::endl;

    for (unsigned int step = 1; step <= n_time_steps; ++step)
      {
        time_step();

        if (step % 10 == 0) // 10ステップ毎に出力
          output_results(timestep_number);

        std::cout << "Time step " << step << " finished, time = " << time << std::endl;
      }
  }
} // namespace ElasticWave2D


int main()
{
  try
    {
      using namespace ElasticWave2D;

      ElasticWave<2> elastic_wave_problem;
      elastic_wave_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "Unknown exception!" << std::endl;
      return 1;
    }

  return 0;
}
