#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <iostream>

using namespace dealii;

int main()
{
  try
  {
    Triangulation<2> triangulation;
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);

    std::ifstream input_file("../simulation1.msh");

    AssertThrow(input_file, ExcMessage("Error: simulation1.msh ファイルを開けませんでした。"));

    grid_in.read_msh(input_file);
    std::cout << "✅ メッシュファイルの読み込みに成功しました。" << std::endl;
    std::cout << "  セル数: " << triangulation.n_active_cells() << std::endl;

    // 追加部分：material_idの表示
    std::cout << "---- 各セルの material_id ----" << std::endl;
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      std::cout << "Cell index: " << cell->index()
                << ", material_id: " << cell->material_id() << std::endl;
    }

  }
  catch (std::exception &exc)
  {
    std::cerr << "❌ 例外発生:\n" << exc.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "❌ 未知の例外が発生しました。" << std::endl;
    return 1;
  }


  return 0;
}
