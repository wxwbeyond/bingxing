#include"mybitset.h"
using namespace std;

int main() 
{
	//��̬�㷨����Ҫ�ֶ����ĳ����е�ȫ�ֳ���
	cout << "�����СΪ" << column_num_c << "����Ԫ�Ӹ���Ϊ" << ek_num_c << "������Ԫ������Ϊ" << et_num_c << endl;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_OMP();  // ��̬λ���洢�ľ���������˹��ȥ
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_OMP: " << (tail - head) * 1000.0 / freq<< " ms" << endl;
}