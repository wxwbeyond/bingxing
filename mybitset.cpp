#include"mybitset.h"

using namespace std;

MyBitSet::MyBitSet()
	{
		head = -1; 
		content = new int[bit_size];
		for (int i = 0; i < bit_size; i++) content[i] = 0;
	}
MyBitSet& MyBitSet:: operator^=(const MyBitSet& b)
	{ 
		for (int i = 0; i < bit_size; i++)
			content[i] ^= b.content[i];
		for (int i = 0; i < bit_size; i++)
		{
			for (int j = 0; j < 32; j++)
			{
				if ((content[i] & (1 << j)))
				{
					head = i * 32 + j;
					return *this;
				}
			}
		}
		head = -1;
		return *this;
	}
MyBitSet& MyBitSet::my_xor_AVX(const MyBitSet& b)
{
		__m256i v_this, v_b;
		int i = 0;
		for (i; i < bit_size - 8; i += 8) 
		{
			v_this = _mm256_loadu_si256((__m256i*) & content[i]);
			v_b = _mm256_loadu_si256((__m256i*) & b.content[i]);
			v_this = _mm256_xor_si256(v_this, v_b);
			_mm256_storeu_si256((__m256i*) & content[i], v_this);
		}
		for (i; i < bit_size; i++)
			content[i] ^= b.content[i];
		for (int i = 0; i < bit_size; i++)
		{
			for (int j = 0; j < 32; j++)
			{
				if ((content[i] & (1 << j)))
				{
					head = i * 32 + j;
					return *this;
				}
			}
		}
		head = -1;
		return *this;
	}
int MyBitSet::test(int index)
	{
		return content[index / 32] & (1 << (index % 32)) ? 1 : 0;  
	}
void MyBitSet::set(int index)
	{ 
		content[index / 32] |= (1 << (index % 32));
	}
bool MyBitSet::any()
	{
	for (int i = 0; i < bit_size; i++)
	{
		if (content[i])
			return true;
	}
		return false;
	}
int find_first_bitset(const bitset<column_num_c>& b)
{
	for (int i = 0; i < column_num_c; i++)
	{
		if (b.test(i))
			return i;
	}
	return -1;
}
void readData_reverse_bitset_c()
{  // 倒序读入数据，读入静态位集
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // 消元子
	ifstream inElimTar(dir + "elimtar.txt");  // 被消元行
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // 用于数据读入
	int lp = -1;
	while (true)  // 读取消元子
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			if (lp == -1)
				lp = column_num_c - inek_loc - 1;
			eks_c[lp].set(column_num_c - inek_loc - 1);
		}
		lp = -1, p_ek++;
		if (inek.empty())
			break;
	}
	while (true)  // 读取被消元行
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			if (lp == -1)
			{
				lp = column_num_c - inet_loc - 1;
				lp_ets_c[p_et] = lp;
			}
			ets_c[p_et].set(column_num_c - inet_loc - 1);
		}
		if (inet.empty())
			break;
		lp = -1;
		p_et++;
	}
	inElimKey.close();
	inElimTar.close();
}
void init_c()
{
	for (int i = 0; i < column_num_c; i++)  // 初始化，处理上一次算法中的残余数据
	{
		eks_c[i] = *(new bitset<column_num_c>);
		lp_eks_c[i] = -1;
	}
	readData_reverse_bitset_c();  // 逆序初始化消元子和被消元行阵列
	cout << "init_complete" << endl;
}
void GrobnerGE_OMP()
{
	int i, j, k;
#pragma omp parallel num_threads(THREAD_NUM), private(i, j, k)
	for (int i = 0; i < column_num_c; i++) // 取每个消元子，对被消元行进行操作，便于并行化
	{
		if (!eks_c[i].test(i)) // 消元子被逆序初始化时满足“行号” = “首项”的条件
		{
#pragma omp barrier
#pragma omp single
			for (size_t j = 0; j < et_num_c; j++)
			{
				if (i == lp_ets_c[j]) // 说明存在对应被消元行
				{
					eks_c[i] = ets_c[j];
					lp_ets_c[j] = -1;
					break;
				}
			}
		}
#pragma omp for schedule(dynamic)
		for (int j = 0; j < et_num_c; j++) // 循环划分并行化
		{
			if (i == lp_ets_c[j]) // 说明存在对应被消元行
			{
				ets_c[j] ^= eks_c[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets_c[j].test(k))
					{
						lp_ets_c[j] = k;
						break;
					}
					if (k == column_num_c - 1)
						lp_ets_c[j] = -1;
				}
			}
		}
	}
}
void readData_reverse_MyB()
{  // 倒序读入数据，读入静态位集
    string inek, inet;
    stringstream ss_inek, ss_inet;
    ifstream inElimKey(dir + "elimkey.txt");  // 消元子
    ifstream inElimTar(dir + "elimtar.txt");  // 被消元行
    int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // 用于数据读入
    int lp = -1;
    while (true)  // 读取消元子
    {
        getline(inElimKey, inek);
        ss_inek = stringstream(inek);
        while (ss_inek >> inek)
        {
            inek_loc = stoi(inek);
            if (lp == -1)
            {
                lp = column_num_c - inek_loc - 1;
                eks[lp].head = lp;
            }
            eks[lp].set(column_num_c - inek_loc - 1);

        }
        lp = -1;  p_ek++;
        if (inek.empty())
            break;
    }
    while (true)  // 读取被消元行
    {
        getline(inElimTar, inet);
        ss_inet = stringstream(inet);
        while (ss_inet >> inet)
        {
            inet_loc = stoi(inet);
            if (lp == -1)
            {
                lp = column_num_c - inet_loc - 1;
                ets[p_et].head = lp;
            }
            ets[p_et].set(column_num_c - inet_loc - 1);
        }
        lp = -1;  p_et++;
        if (inet.empty()) break;
    }
    inElimKey.close();
    inElimTar.close();
}
void init_MyB()
{
	readData_reverse_MyB();  // 逆序初始化消元子和被消元行阵列
	cout << "init_complete" << endl;
}
