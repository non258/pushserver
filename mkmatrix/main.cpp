#include <fstream>
#include <iostream>
#include <vector>

int main()
{
  const int size = 10000000;
  const char *fileNameval = "../7_25val.txt";
  const char *fileNamerow = "../7_25row.txt";
  const char *fileNamecol = "../7_25col.txt";

  int rowcount = 0;
  int colcount = 0;

  std::ofstream ofsval(fileNameval);
  std::ofstream ofsrow(fileNamerow);
  std::ofstream ofscol(fileNamecol);
  if (!ofsval || !ofsrow || !ofscol)
  {
    std::cout << "ファイルが開けませんでした。" << std::endl;
    std::cin.get();
    return 0;
  }

  for (int i = 0; i < size; i++)
  {
    std::cout << size-i << std::endl;
    if (i == 0)
    {
      ofsval << "25" << std::endl;
      ofsval << "-1" << std::endl;

      ofsrow << "0" << std::endl;
      ofsrow << "25" << std::endl;
      rowcount += 2;

      ofscol << "0" << std::endl;
      ofscol << "1" << std::endl;
    }
    else if (i == size-1)
    {
      ofsval << "-1" << std::endl;
      ofsval << "25" << std::endl;

      std::string str = std::to_string(rowcount+2);
      ofsrow << str << std::endl;
      rowcount += 2;

      ofscol << std::to_string(colcount) << std::endl;
      ofscol << std::to_string(colcount+1) << std::endl;
      colcount++;
    }
    else
    {
      ofsval << "-1" << std::endl;
      ofsval << "25" << std::endl;
      ofsval << "-1" << std::endl;

      std::string str = std::to_string(rowcount+3);
      ofsrow << str << std::endl;
      rowcount += 3;

      ofscol << std::to_string(colcount) << std::endl;
      ofscol << std::to_string(colcount+1) << std::endl;
      ofscol << std::to_string(colcount+2) << std::endl;
      colcount++;
    }
  }
}
