// deal.IIバージョン: 9.5以降推奨（FE_SimplexPと三角形メッシュ対応）
// 三角形要素に対応した2D弾性波動方程式の陽的時間積分
//P波のみ　or　S波のみで出力してそうなっているかをparaviewで確認するよう 
//const bool excite_p_wave = false;  // ← falseにすればS波のみ励起　ここで調整　
//これで分離できると思ったが駄目だった。おそらく初期条件をある頂点での変位としているためだろう

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
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_q_generic.h>

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
    // Vp, Vs, rho の物性値を material_id に応じて用意（例）

    static constexpr double Vp_ice     = 3850.0 * 1000.0;   // [mm/s]
    static constexpr double Vs_ice     = 1950.0 * 1000.0;   // [mm/s]
    static constexpr double rho_ice    = 917.0 / 1.0e9;     // [kg/mm^3]

    static constexpr double Vp_shale   = 3500.0 * 1000.0;
    static constexpr double Vs_shale   = 1800.0 * 1000.0;
    static constexpr double rho_shale  = 2400.0 / 1.0e9;

    static constexpr double Vp_gneiss  = 5500.0 * 1000.0;
    static constexpr double Vs_gneiss  = 3000.0 * 1000.0;
    static constexpr double rho_gneiss = 2650.0 / 1.0e9;


    static double rho(const types::material_id id)
    {
        switch(id)
        {
        case 19: return rho_ice;
        case 20: return rho_shale;
        case 21: return rho_gneiss;
        default: return 1.0; // default値
        }
    }

    static double lambda(const types::material_id id)
    {
        // λ = ρ(Vp² - 2Vs²)
        const double rho_val = rho(id);
        const double Vp = Vp_of(id);
        const double Vs = Vs_of(id);
        return rho_val * (Vp*Vp - 2.0*Vs*Vs);
    }

    static double mu(const types::material_id id)
    {
        // μ = ρ * Vs²
        return rho(id) * Vs_of(id)*Vs_of(id);
    }

    private:
    static double Vp_of(const types::material_id id)
    {
        switch(id)
        {
        case 19: return Vp_ice;
        case 20: return Vp_shale;
        case 21: return Vp_gneiss;
        default: return 1.0;
        }
    }
    static double Vs_of(const types::material_id id)
    {
        switch(id)
        {
        case 19: return Vs_ice;
        case 20: return Vs_shale;
        case 21: return Vs_gneiss;
        default: return 1.0;
        }
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
    , fe(FE_SimplexP<dim>(1), dim) // 三角形要素に変更
    , time(0.)
    , time_step_size(1e-4)
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
    //以下のvectorToolsの行を省くとノイマン条件になり省かないとディリクレ条件になる
    // VectorTools::interpolate_boundary_values(dof_handler, 22, Functions::ZeroFunction<dim>(dim), constraints);
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
    QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points = quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {

      // std::cout << "Cell " << cell->index() << " material_id = " << cell->material_id() << std::endl;
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
    QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
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

        const Point<dim> source_center(280000.0, 100000.0);

        double min_distance = std::numeric_limits<double>::max();
        types::global_dof_index nearest_dof = numbers::invalid_dof_index;
        Point<dim> nearest_point;

        // 最も近い頂点を探す
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            const unsigned int n_vertices = cell->n_vertices();

            for (unsigned int vertex_index = 0; vertex_index < n_vertices; ++vertex_index)
            {
                const Point<dim> &p = cell->vertex(vertex_index);
                const double dist = p.distance(source_center);
                if (dist < min_distance)
                {
                    nearest_dof = cell->vertex_dof_index(vertex_index, 0); // 任意の成分
                    min_distance = dist;
                    nearest_point = p;
                }
            }
        }

        Assert(nearest_dof != numbers::invalid_dof_index, ExcMessage("No closest vertex found!"));

        // ======== 選択：P波 or S波 ========
        const bool excite_p_wave = false;  // ← falseにすればS波のみ励起

        // 最近接点からのベクトル方向（単位ベクトル）
        Tensor<1, dim> direction = nearest_point - source_center;
        direction /= direction.norm();

        if (!excite_p_wave)
        {
            // S波のため、伝搬方向に垂直なベクトルを使う（2D限定の単純な回転）
            std::swap(direction[0], direction[1]);
            direction[0] *= -1; // 90度回転（時計回り）
        }

        // 最近接点の自由度に方向付きの変位を与える（ベクトル成分ごとに）
        for (unsigned int d = 0; d < dim; ++d)
        {
            
            for (const auto &cell : dof_handler.active_cell_iterators())
            {
                for (unsigned int v = 0; v < cell->n_vertices(); ++v)
                {
                    if (cell->vertex(v).distance(nearest_point) < 1e-12)
                    {
                        const types::global_dof_index dof = cell->vertex_dof_index(v, d);
                        solution_n(dof) = 0.1 * direction[d];
                    }
                }
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
    {
        double term = mass_matrix_diagonal_inverse[i] * system_rhs[i] * time_step_size * time_step_size
                    + 2.0 * solution_n[i] - solution_nm1[i];

        solution_np1[i] = term;
    }

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

        // displacement ベクトル場をベクトルとして追加（重要）
        data_out.add_data_vector(solution_n, "displacement",
            DataOut<dim>::type_dof_data,
            std::vector<DataComponentInterpretation::DataComponentInterpretation>(
                dim, DataComponentInterpretation::component_is_part_of_vector));

        // Q点設定
        const QGaussSimplex<dim> quadrature(fe.degree + 1);
        FEValues<dim> fe_values(fe, quadrature,
                                update_gradients | update_quadrature_points);

        const unsigned int n_q_points = quadrature.size();
        Vector<double> div_values(dof_handler.n_dofs());
        Vector<double> curl_values(dof_handler.n_dofs());

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            fe_values.reinit(cell);
            const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
            std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);

            Vector<double> local_solution(dofs_per_cell);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                local_solution[i] = solution_n(local_dof_indices[i]);

            Tensor<2, dim> avg_grad;
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                Tensor<2, dim> grad_u;
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                    const unsigned int comp = fe.system_to_component_index(i).first;
                    grad_u[comp] += local_solution[i] * fe_values.shape_grad(i, q);
                }
                avg_grad += grad_u;
            }
            avg_grad /= n_q_points;

            // 各DOFに div と curl を代入（全DOFに同じ値を代入）
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const auto gdof = local_dof_indices[i];
                double div = trace(avg_grad); // = grad_00 + grad_11
                div_values[gdof] = div;

                if constexpr (dim == 2)
                {
                    double curl = avg_grad[1][0] - avg_grad[0][1];
                    curl_values[gdof] = curl;
                }
            }
        }

        data_out.add_data_vector(div_values, "div_u");
        if constexpr (dim == 2)
            data_out.add_data_vector(curl_values, "curl_u");

        data_out.build_patches();

        const std::string filename = "solution-" + std::to_string(timestep) + ".vtu";
        std::ofstream output(filename);
        data_out.write_vtu(output);

        std::cout << "VTU file written: " << filename << std::endl;
    }







    template <int dim>
    void ElasticWave<dim>::run()
    {
    // GMSHメッシュを読み込む部分
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);

    std::ifstream input_file("../simulation1.msh");
    Assert(input_file, ExcFileNotOpen("simulation1.msh"));

    grid_in.read_msh(input_file);


/////////////////////////////////////////////

// 🔽 境界IDの確認
    {
    std::set<types::boundary_id> boundary_ids;
    for (const auto &cell : triangulation.active_cell_iterators())
        for (unsigned int f = 0; f < cell->n_faces(); ++f)
        if (cell->face(f)->at_boundary())
            boundary_ids.insert(cell->face(f)->boundary_id());

    std::cout << "Found boundary IDs in mesh: ";
    for (const auto &id : boundary_ids)
        std::cout << static_cast<unsigned int>(id) << " ";
    std::cout << std::endl;
    }

///////////////////////////////////////////////



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