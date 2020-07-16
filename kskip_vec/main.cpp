#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

//std::vector<double> mvm(std::vector<std::vector<double>> A, std::vector<double> b, int size)
//{
//  std::vector<double> result(size, 0);
//  for (int i = 0; i < size; i++)
//    for (int j = 0; j < size; j++)
//      result[i] += A[i][j] * b[j];
//
//  return result;
//}
//
template <typename I, std::enable_if_t<std::is_integral_v<I> && std::is_arithmetic_v<I>, std::nullptr_t> = nullptr, typename T, std::enable_if_t<std::is_arithmetic_v<T>, std::nullptr_t> = nullptr>
void mvm(T *val, I *col, I *row, T *b, I size, T *Ax)
{
  // std::cout << "----------end----------" << std::endl;
  // for (int i = 0; i < size+1; i++)
  // {
  //   std::cout << row[i] << std::endl;
  // }
  for (int i = 0; i < size; i++)
  {
    for (int j = row[i]; j < row[i+1]; j++)
    {
      // std::cout << row[i] << ":" << row[i+1] << std::endl;
      Ax[i] += val[j] * b[col[j]];
      // std::cout << result[i] << std::endl;
    }
  }
  // std::cout << "----------end----------" << std::endl;
}

template <typename I, std::enable_if_t<std::is_integral_v<I> && std::is_arithmetic_v<I>, std::nullptr_t> = nullptr, typename T, std::enable_if_t<std::is_arithmetic_v<T>, std::nullptr_t> = nullptr>
void mvm(T *val, I *col, I *row, T *b, I size, T *Ax, I s)
{
  // std::cout << "----------end----------" << std::endl;
  // for (int i = 0; i < size+1; i++)
  // {
  //   std::cout << row[i] << std::endl;
  // }
  for (int i = 0; i < size; i++)
  {
    for (int j = row[i]; j < row[i+1]; j++)
    {
      // std::cout << row[i] << ":" << row[i+1] << std::endl;
      Ax[i] += val[j] * b[s*size+col[j]];
      // std::cout << result[i] << std::endl;
    }
  }
  // std::cout << "----------end----------" << std::endl;
}
//
double ip(double *a, int sa, double *b, int sb, int size)
{
  double result = 0;
  for (int i = 0; i < size; i++)
  {
    result += a[sa*size+i] * b[sb*size+i];
//    std::cout << result << std::endl;
  }
  return result;
}
//
//std::vector<std::vector<double>> matrixpow (std::vector<std::vector<double>> A, int n, int size)
//{
//  for (int i = 0; i < size; i ++)
//  {
//    for (int j = 0; j < size; j++)
//    {
//      A[i][j] = pow(A[i][j], n);
//    }
//  }
//
//  return A;
//}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, std::nullptr_t> = nullptr>
void inputdata(std::string filename, T *arr, int size)
{
//  std::cout << "ok" << std::endl;
  std::ifstream ifs(filename);
  std::string str;

  if (ifs.fail())
    std::cout<< "Failded to open vec file." << std::endl;

  int count = 0;
  while (getline(ifs, str))
  {
    arr[count] = std::stod(str);
    count++;
    if (size <= count)
      break;
  }
}

//std::vector<double> inputvec(std::string filename, int size)
//{
//  std::vector<double> input(size, 0);
//  std::ifstream ifs(filename);
//  std::string str;
//
//  if (ifs.fail())
//    std::cout<< "Failded to open vec file." << std::endl;
//
//  int count = 0;
//  while (getline(ifs, str))
//  {
//    input[count] = std::stod(str);
//    count++;
//    if (size <= count)
//      break;
//  }
//
//  return input;
//}
//
//std::vector<int> inputvecint(std::string filename, int size)
//{
//  std::vector<int> input(size, 0);
//  std::ifstream ifs(filename);
//  std::string str;
//
//  if (ifs.fail())
//    std::cout<< "Failded to open vec file." << std::endl;
//
//  int count = 0;
//  while (getline(ifs, str))
//  {
//    input[count] = std::stod(str);
//    count++;
//    if (size <= count)
//      break;
//  }
//
//  return input;
//}
//
//std::vector<std::vector<double>> inputmatrix(std::string filename, int size)
//{
//  std::vector<std::vector<double>> input(size, std::vector<double>(size, 0));
//  std::ifstream ifs(filename);
//  std::string str;
//
//  if (ifs.fail())
//    std::cout << "Failded to open matrix file." << std::endl;
//
//  int count1 = 0;
//  int count2 = 0;
//  while (getline(ifs, str))
//  {
//    input[count1][count2] = std::stod(str);
//    count2++;
//    if (size <= count2)
//    {
//      count2 = 0;
//      count1++;
//      if (size <= count1)
//        break;
//    }
//  }
//
//  return input;
//}
//
double norm(double *vec, int size)
{
  double s = 0.0;

  for (int i = 0; i < size; i++)
  {
    s += vec[i] * vec[i];
  }

  return sqrt(s);
}

// void printvector(std::vector<double> vec, int size)
// {
//   for (int i = 0; i < size; i++)
//     std::cout << vec[i] << " " << std::endl;
// }
//
// void printmatrix(std::vector<std::vector<double>> matrix, int size)
// {
//   for (int i = 0; i < size; i++)
//   {
//     for (int j = 0; j < size; j++)
//     {
//       std::cout << matrix[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }
// }

int main() {
  int size = 1000000;
  int val_size = 2999998;
  std::string filename = "6_20005";
  int k = 4;
  double epsilon = 1*std::pow(10,-10);
  double d = 0;
  int step = 0;

  // std::vector<std::vector<double>> A = inputmatrix("../3a.txt", size);

  double *val = new double[val_size]();
  int *col = new int[val_size]();
  int *row = new int[size+1]();
  double *b = new double[size]();
  double *x = new double[size]();

inputdata("../data/" + filename + "val.txt", val, val_size);
//  for (int i = 0; i < val_size; i++)
//  {
//    std::cout << val[i] << std::endl;
//  }
//  return 0;
  inputdata("../data/" + filename + "col.txt", col, val_size);
//  for (int i = 0; i < val_size; i++)
//  {
//    std::cout << col[i] << std::endl;
//  }
//  return 0;
  inputdata("../data/" + filename + "row.txt", row, size+1);
//  for (int i = 0; i < size+1; i++)
//  {
//    std::cout << row[i] << std::endl;
//  }
//  return 0;
  inputdata("../data/" + std::to_string(size) + "b.txt", b, size);
//  for (int i = 0; i < size; i++)
//  {
//    std::cout << b[i] << std::endl;
//  }
//  return 0;
//
//
  std::chrono::steady_clock::time_point start;
  double time;
//
  double *Ax = new double[size]();
  mvm(val, col, row, x, size, Ax);
//  for (int i = 0; i < size; i++)
//  {
//    std::cout << Ax[i] << std::endl;
//  }
//  return 0;

//  // printvector(Ax, size);
//
//  //ok
//  std::vector<std::vector<double>> Ar(k+2, std::vector<double>(size, 0));
  double *Ar = new double[(k+2)*size]();
//  std::cout << (k+2) * size << std::endl;
//  for (int i = 0; i < (k+2)*size; i++)
//  {
//    std::cout << Ar[i] << std::endl;
//  }
//  return 0;
//  // r = b-Ax
  for (int i = 0; i < size; i++)
    Ar[i] = b[i] - Ax[i];
//  for (int i = 0; i < size; i++)
//  {
//    std::cout << Ar[i] << std::endl;
//  }
//  return 0;
//
//  //ok
//  std::vector<std::vector<double>> Ap(k+3, std::vector<double>(size, 0));
  double *Ap = new double[(k+2)*size]();
//  // p = r
  for (int i = 0; i < size; i++)
    Ap[i] = Ar[i];
//  for (int i = 0; i < size; i++)
//  {
//    std::cout << Ap[i] << std::endl;
//  }
//  return 0;
//
//  // printmatrix(Ap, size);
//
//  std::vector<double>a(2*k+2, 0);
//  std::vector<double>f(2*k+4, 0);
//  std::vector<double>c(2*k+2, 0);
  double *a = new double[2*k+2]();
  double *f = new double[2*k+4]();
  double *c = new double[2*k+2]();
//
//  std::vector<double>residual(size+1, 0);
  double *residual = new double[size+1]();
//
  for (int i = 0; i < size; i++)
  {
    std::cout << "step: " << i+1 << std::endl;
    //ok
    residual[i] = norm(Ar, size) / norm(b, size);
//    std::cout << residual[i] << std::endl;
//    return 0;
    // std::cout << residual[i] << std::endl;
    // return 0;
    // std::cout << residual[i] << " : " << epsilon << std::endl;
    if (residual[i] < epsilon)
    {
      step = i+1;
      break;
    }

    // ok
    for (int j = 1; j < k+1; j++)
    {
//      std::vector<double> tmpAr = mvmA(val, col, row, size, Ar[j-1]);
      double *tmpAr = new double[size]();
//      void mvm(T *val, I *col, I *row, T *b, I size, T *Ax, I s)
      mvm(val, col, row, Ar, size, tmpAr, j-1);
//      for (int l = 0; l < size; l++)
//        std::cout << tmpAr[l] << std::endl;
//      return 0;
      for (int l = 0; l < size; l++)
      {
        Ar[j*size+l] = tmpAr[l];
//        std::cout << Ar[j*size+l] << std::endl;
      }
    }
//    return 0;

    // printmatrix(Ar, size);
    // return 0;

    //ok
    for (int j = 1; j < k+2; j++)
    {
//      std::vector<double> tmpAp = mvmA(val, col, row, size, Ap[j-1]);
      double *tmpAp = new double[size]();
      mvm(val, col, row, Ap, size, tmpAp, j-1);
//      for (int l = 0; l < size; l++)
//        std::cout << tmpAp[l] << std::endl;
//      return 0;

      for (int l = 0; l < size; l++)
        Ap[j*size+l] = tmpAp[l];
    }
    // for (int j = 0; j < k+3; j++)
    // {
    //   for (int l = 0; l < size; l++)
    //   {
    //     std::cout << Ap[j][l] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // return 0;

//    for (int j = 0; j < 2*k+1; j++)
//    {
//      for (int l = 0; l < size; l++)
//      {
//        std::cout << Ar[j*size+l] << " ";
//      }
//      std::cout << std::endl;
//    }
//    return 0;

    // ok
    for (int j = 0; j < 2*k+1; j += 2)
    {
      int jj = j/2;
      a[j] = ip(Ar, jj, Ar, jj, size);
//      std::cout << a[j] << std::endl;
      a[j+1] = ip(Ar, jj, Ar, jj+1, size);
//      std::cout << a[j+1] << std::endl;
//      for (int l = 0; l < size; l++)
//        std::cout << Ar[(jj+1)*size+l] << " ";
//      std::cout << std::endl;
    }
//    return 0;
    // for (int j = 0; j < a.size(); j++)
    //   std::cout << a[j] << " ";
    // std::cout << std::endl;
    // return 0;

    // ok
    for (int j = 0; j < 2*k+3; j += 2)
    {
      int jj = j/2;
      f[j] = ip(Ap, jj, Ap, jj, size);
      f[j+1] = ip(Ap, jj, Ap, jj+1, size);
    }
    // for (int j = 0; j < f.size(); j++)
    //   std::cout << f[j] << " ";
    // std::cout << std::endl;
    // return 0;

    // ok
    for (int j = 0; j < 2*k+1; j += 2)
    {
      int jj = j/2;
      c[j] = ip(Ar, jj, Ap, jj, size);
      c[j+1] = ip(Ar, jj, Ap, jj+1, size);
    }
    // for (int j = 0; j < c.size(); j++)
    //   std::cout << c[j] << " ";
    // std::cout << std::endl;
    // return 0;

    // ok
//    alpha = a[0] / f[1]
    double alpha = a[0] / f[1];
    // std::cout << alpha << std::endl;
    // return 0;
    //ok
//    beta = alpha**2 * f[2] / a[0] - 1
    double beta = pow(alpha,2) * f[2] / a[0] - 1;
    // std::cout << beta << std::endl;
    // return 0;

    // ok
//    self.x += alpha*Ap[0]
    for (int j = 0; j < size; j++)
      x[j] += alpha * Ap[j];
    //  for (int j = 0; j < x.size(); j++)
    //    std::cout << x[j] << " ";
    // std::cout << std::endl;
    // return 0;

    // ok
//    Ar[0] -= alpha*Ap[1]
    for (int j = 0; j < size; j++)
      Ar[j] -= alpha * Ap[1*size+j];
    // for (int j = 0; j < Ar.size(); j++)
    // {
    //   for (int l = 0; l < Ar[j].size(); l++)
    //   {
    //     std::cout << Ar[j][l] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // return 0;

    // ok
//    Ap[0] = Ar[0] + beta*Ap[0]
    for (int j = 0; j < size; j++)
      Ap[j] = Ar[j] + beta * Ap[j];
    // for (int j = 0; j < Ap.size(); j++)
    // {
    //   for (int l = 0; l < Ap[j].size(); l++)
    //   {
    //     std::cout << Ap[j][l] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // return 0;


    // ok
//    Ap[1] = self.A.dot(Ap[0])
    {
//      std::vector<double> tmpAp = mvmA(val, col, row, size, Ap[0]);
      double *tmpAp = new double[size]();
      mvm(val, col, row, Ap, size, tmpAp);
      for (int j = 0; j < size; j++)
        Ap[1*size+j] = tmpAp[j];
    }
    // for (int j = 0; j < Ap.size(); j++)
    // {
    //   for (int l = 0; l < Ap[j].size(); l++)
    //   {
    //     std::cout << Ap[j][l] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // return 0;


//    for j in range(k):
    for (int j = 0; j < k; j++)
    {
//    for l in range(0, 2*(k-j)+1):
      for (int l = 0; l < 2*(k-j)+1; l++)
      {
        //std::cout << "out" << 2*(k-j)+1 << std::endl;
        //ok
        a[l] += alpha*(alpha*f[l+2] - 2*c[l+1]);
        // std::cout << "l:" << l << std::endl;
        // std::cout << a[l] << std::endl;
        d = c[l] - alpha*f[l+1];
        // std::cout << d << std::endl;
        c[l] = a[l] + d*beta;
        // std::cout << c[l] << std::endl;
        f[l] = c[l] + beta*(d + beta*f[l]);
        // std::cout << f[l] << std::endl;
      }
      alpha = a[0] / f[1];
      // std::cout << alpha << std::endl;
      beta = pow(alpha, 2) * f[2] / a[0] - 1;
//        x += alpha*Ap[0]
      for (int m = 0; m < size; m++)
        x[m] += alpha*Ap[m];
//        Ar[0] -= alpha*Ap[1]
      for (int m = 0; m < size; m++)
        Ar[m] -= alpha*Ap[1*size+m];
//        Ap[0] = Ar[0] + beta*Ap[0]
      for (int m = 0; m < size; m++)
        Ap[m] = Ar[m] + beta*Ap[m];
//        Ap[1] = self.A.dot(Ap[0])
      {
//        std::vector<double> tmpAp = mvmA(val, col, row, size, Ap[0]);
        double *tmpAp = new double[size]();
        mvm(val, col, row, Ap, size, tmpAp);
        for (int m = 0; m < size; m++)
          Ap[1*size+m] = tmpAp[m];
      }
    }
  }

  std::chrono::steady_clock::time_point const end = std::chrono::steady_clock::now();
  time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  // std::printf("time : %g [ms]\n", time);

  // std::cout << "-----residual-----" << std::endl;
  // for (int i = 0; i < size; i++)
  //   std::cout << residual[i] << std::endl;
  // std::cout << "-----result-----" << std::endl;
  // for (int i = 0; i < size; i++)
  //   std::cout << x[i] << std::endl;

  const std::string fileName = "../data/" + filename + "k" + std::to_string(k) + "_result.txt";

  std::ofstream ofs(fileName);
  if (!ofs)
  {
    std::cout << "ファイルが開けませんでした。" << std::endl;
    std::cin.get();
    return 0;
  }

  ofs << "-----step-----" << std::endl;
  ofs << step << std::endl;
  ofs << "-----time-----" << std::endl;
  ofs << time << "[ms]" << std::endl;
  ofs << "-----residual-----" << std::endl;
  for (int i = 0; i < step; i++)
    ofs << residual[i] << std::endl;
  // std::cout << "-----result-----" << std::endl;
  //  for (int i = 0; i < size; i++)
  //    std::cout << x[i] << std::endl;
}
