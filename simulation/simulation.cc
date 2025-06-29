#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim>
class WaveEquation
{
public:
  WaveEquation();
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void apply_initial_conditions();
  void time_step();
  void output_results(const unsigned int timestep) const;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  Vector<double> solution;        // 変位 u^n
  Vector<double> solution_old;    // 変位 u^{n-1}
  Vector<double> solution_new;    // 変位 u^{n+1}
  Vector<double> system_rhs;

  // 物性パラメータ：P波速度, S波速度, 密度を領域ごとに保存
  std::map<unsigned int, double> p_wave_speed;
  std::map<unsigned int, double> s_wave_speed;
  std::map<unsigned int, double> density;

  double dt;  // 時間刻み
  unsigned int n_time_steps;

  // ここではシンプルな陽的時間積分
};

template <int dim>
WaveEquation<dim>::WaveEquation()
  : fe(1)
  , dof_handler(triangulation)
  , dt(0.001)
  , n_time_steps(1000)
{}

template <int dim>
void WaveEquation<dim>::make_grid()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("simulation1.msh");
  Assert(dim == 2, ExcNotImplemented());
  grid_in.read_msh(input_file);

  std::cout << "Cells: " << triangulation.n_cells() << std::endl;
  std::cout << "Vertices: " << triangulation.n_vertices() << std::endl;

  triangulation.refine_global(0); // ここは必要に応じて調整

  std::cout << "Grid loaded." << std::endl;
}

template <int dim>
void WaveEquation<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  solution_old.reinit(dof_handler.n_dofs());
  solution_new.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  // 物性パラメータ設定（material_id = 0,1,2 に対応）
  // P波速度[m/s], S波速度[m/s], 密度[kg/m3] は適宜実測値に置き換えてください
  p_wave_speed[0] = 3800.0;   // 氷層
  s_wave_speed[0] = 1900.0;
  density[0]      = 900.0;

  p_wave_speed[1] = 5000.0;   // 両脇層
  s_wave_speed[1] = 3000.0;
  density[1]      = 2700.0;

  p_wave_speed[2] = 6000.0;   // 下部gneiss層
  s_wave_speed[2] = 3500.0;
  density[2]      = 2700.0;
}

template <int dim>
void WaveEquation<dim>::apply_initial_conditions()
{
  solution = 0;
  solution_old = 0;
  solution_new = 0;
}

template <int dim>
void WaveEquation<dim>::assemble_system()
{
  // ここでは簡易的に system_rhs に波源の力を設定する例を示す

  // 波源位置の範囲（x座標, y座標）を指定
  // 氷層はx=0..500m, y=0..240m
  // 両脇層はx=-30..0と500..530m, y=0..240m
  // 下部層はx=0..500m, y=-30..0m
  //
  // 波源は氷層の右側の層の右上端から右に水平に100m部分、例えば (530..630, 210..240) に置く想定
  // ここでは簡単に (530 <= x <= 630, 210 <= y <= 240) に非ゼロの力を与える

  // まずrhsを0に
  system_rhs = 0;

  // 各自由度の位置を取得
  std::vector<Point<dim>> support_points(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, support_points);

  // 力の強さ
  const double force_magnitude = 1e6; // 適宜調整

  for (unsigned int i=0; i < dof_handler.n_dofs(); ++i)
  {
    // 変位はdim成分なので、iがどの成分か判定する必要あり
    const unsigned int component_i = fe.system_to_component_index(i).first;
    // 水平方向成分（component 0）にのみ力を加える想定
    if (component_i == 0)
    {
      const Point<dim> &p = support_points[i];
      if (p[0] >= 530 && p[0] <= 630 &&
          p[1] >= 210 && p[1] <= 240)
      {
        system_rhs[i] = force_magnitude;
      }
    }
  }
}

template <int dim>
void WaveEquation<dim>::time_step()
{
  // ここでは単純な陽的中央差分法の例（非常に簡易、安定条件注意）

  // 解の自由度数
  const unsigned int n_dofs = dof_handler.n_dofs();

  // 運動方程式の簡易離散化
  // M (u^{n+1} - 2 u^{n} + u^{n-1}) / dt^2 = F - K u^n
  // ここでは質量行列Mと剛性行列Kを仮定せず、単純化のため
  // system_rhsに外力を代入、M,Kを単位行列として近似

  for (unsigned int i=0; i<n_dofs; ++i)
  {
    solution_new[i] = 2*solution[i] - solution_old[i] + dt*dt*system_rhs[i];
  }

  solution_old = solution;
  solution = solution_new;
}

template <int dim>
void WaveEquation<dim>::output_results(const unsigned int timestep_number) const
{
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    // displacement vectorを出力（ベクトル値）
    data_out.add_data_vector(solution, displacement_names, DataOut<dim>::type_dof_data);

    // material_idはセルデータなので、unsigned int → double のベクトルに変換して渡す
    std::vector<double> material_ids_double(triangulation.n_active_cells());
    unsigned int idx = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
    {
        material_ids_double[idx] = static_cast<double>(cell->material_id());
        ++idx;
    }
    // セルデータとして追加するために第3引数でtype_cell_dataを指定
    data_out.add_data_vector(material_ids_double, "material_id", DataOut<dim>::type_cell_data);

    data_out.build_patches();

    // ファイル名を作成（例: solution-000.vtu）
    std::string filename = "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);

    std::cout << "Output results to " << filename << std::endl;
}


template <int dim>
void WaveEquation<dim>::run()
{
  make_grid();
  setup_system();
  apply_initial_conditions();

  for (unsigned int timestep=0; timestep<n_time_steps; ++timestep)
  {
    assemble_system();
    time_step();

    if (timestep % 10 == 0)
      output_results(timestep);
  }
}

int main()
{
  try
  {
    WaveEquation<2> wave_equation_2d;
    wave_equation_2d.run();
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
