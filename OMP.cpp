#include"mybitset.h"
using namespace std;

int main() 
{
	//静态算法则需要手动更改程序中的全局常量
	cout << "矩阵大小为" << column_num_c << "，消元子个数为" << ek_num_c << "，被消元行行数为" << et_num_c << endl;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_OMP();  // 动态位集存储的矩阵的特殊高斯消去
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_OMP: " << (tail - head) * 1000.0 / freq<< " ms" << endl;
}