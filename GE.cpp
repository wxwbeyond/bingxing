#include"mybitset.h"
using namespace std;
using namespace boost;


//��������
dynamic_bitset<> get_head(size_t index);
void GrobnerGE();
void db_GrobnerGE();  
void readData_bitset();
void readData_reverse_bitset();



vector<dynamic_bitset<>> eks;
vector<int> lp_eks;
dynamic_bitset<>* ets = nullptr;
int* lp_ets = nullptr;

//�߳����ݽṹ
typedef struct
{
	int t_id;  // �߳� id
	int tasknum;  // ��������
}PT_EliminationParam;
typedef struct
{
	int t_id; //�߳� id
	size_t index;
}PT_GetHeadParam;
//�ź�������
sem_t sem_gethead;
sem_t sem_workerstart[THREAD_NUM]; // ÿ���߳����Լ�ר�����ź���
sem_t sem_workerend[THREAD_NUM];

sem_t sem_leader;
sem_t sem_Divsion[THREAD_NUM - 1];
sem_t sem_Elimination[THREAD_NUM - 1];

//barrier����
pthread_barrier_t barrier_gethead;
pthread_barrier_t _Barrier_Elimination;
pthread_barrier_t _Barrier_0_Operation;

//���������ͻ���������
pthread_mutex_t _mutex;
pthread_cond_t _cond;

//�̺߳���

void* PT_Elimination(void* param)
{
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	bool find_flag = false;
	dynamic_bitset<>* temp = nullptr;
	for (int i = t_id; i < et_num; i += THREAD_NUM)
	{
		while (ets_c[i].any())
		{
			find_flag = false;
			size_t index = ets[i].find_first();
			for (int j = 0; j < eks.size(); j++)  // R[lp(E[i])
			{
				if (index == eks[j].find_first())
				{
					ets[i] ^= eks[j];
					find_flag = true;
					break;
				}
			}
			if (!find_flag)
			{
				pthread_mutex_lock(&_mutex);  // ʹ�û�������֤���ݼ�ͬ��
				eks.push_back(ets[i]);
				pthread_mutex_unlock(&_mutex);  // ����
				break;
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}
unordered_set<int> ek_set;
void* PT_Elimination_Static(void* param)
{
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks_c[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks_c[i] = ets_c[j];
						lp_ets_c[j] = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id; j < et_num_c; j += THREAD_NUM)  // ѭ�����ֲ��л�
		{
			if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
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
	pthread_exit(nullptr);
	return nullptr;
}
void* PT_Elimination_Static_Block(void* param)
{
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	int tasknum = t_id == THREAD_NUM - 1 ? et_num_c % (THREAD_NUM - 1) : et_num_c / (THREAD_NUM - 1);
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks_c[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks_c[i] = ets_c[j];
						lp_ets_c[j] = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id * tasknum, k = 0; k < tasknum; j++, k++)  // �黮�ֲ��л�
		{
			if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets_c[j] ^= eks_c[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets_c[j].test(k))
					{
						lp_ets_c[j] = k;
						break;
					}
					if (k == column_num_c - 1) lp_ets_c[j] = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}
void* PT_Elimination_Static_AVX(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id; j < et_num_c; j += THREAD_NUM)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}
void* PT_Elimination_Static_Block_AVX(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	int tasknum = t_id == THREAD_NUM - 1 ? et_num_c % (THREAD_NUM - 1) : et_num_c / (THREAD_NUM - 1);
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id * tasknum, k = 0; k < tasknum; j++, k++)  // �黮�ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}
void* watch()
{
	pthread_mutex_lock(&_mutex);
	while (true)
		pthread_cond_wait(&_cond, &_mutex);
	pthread_mutex_unlock(&_mutex);
	return nullptr;
}
//PThread�㷨
void PT_GrobnerGE()
{
	pthread_cond_init(&_cond, nullptr);
	pthread_mutex_init(&_mutex, nullptr);
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		pthread_join(handles[t_id], nullptr);
	}
	pthread_mutex_destroy(&_mutex);
	pthread_cond_destroy(&_cond);
	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}
void PT_GrobnerGE_Static()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static, (void*)&param[t_id]);
	}
	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}
void PT_GrobnerGE_Static_Block()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_Block, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}
void PT_GrobnerGE_Static_AVX()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_AVX, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}
void PT_GrobnerGE_Static_Block_AVX()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_Block_AVX, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}

//���㸨������
dynamic_bitset<> get_head(size_t index)
{  // ��ȡ����
	for (int i = 0; i < eks.size(); i++)  // R[lp(E[i])
	{
		if (eks[i].find_first() == index) 
			return eks[i];
	}
	return dynamic_bitset<>(0);
}
dynamic_bitset<> get_head_lp(size_t index)
{  // ��ȡ����
	for (int i = 0; i < eks.size(); i++)  // R[lp(E[i])
	{
		if (index == lp_eks[i])
			return eks[i];
	}
	return dynamic_bitset<>(0);
}
int get_head_lp_c(size_t index)
{  // ��ȡ����
	for (int i = 0; i < column_num_c; i++)
	{
		cout << index << " " << i << endl;
		if (index == lp_eks_c[i])
			return i;
		else if (lp_eks_c[i] == -1)
			break;
	}
	return -1;
}
void xor_static_AVX(bitset<column_num_c>& b1, const bitset<column_num_c>& b2)
{  // b1 ^= b2����֤b1��b2�ȳ�
	for (size_t _Wpos = 0; _Wpos < column_num_c; ++_Wpos)
	{
		b1.set(_Wpos, b1._Subscript(_Wpos) ^ b2._Subscript(_Wpos));
	}
}
//���㸨������-PThread
void* PT_Rotation_gethead_Inthread(void* param) {
	PT_GetHeadParam* pa = (PT_GetHeadParam*)param;
	int t_id = pa->t_id;
	size_t index = pa->index;
	dynamic_bitset<>* ret = nullptr;
	for (int i = t_id; i < eks.size(); i += THREAD_NUM)  // R[lp(E[i])
	{
		if (index == eks[i].find_first())
		{
			ret = new dynamic_bitset<>(eks[i]);
			break;
		}
	}
	sem_post(&sem_gethead);// ֪ͨ leader, �������ȥ����
	pthread_exit(ret);
	return ret;
}
dynamic_bitset<> PT_Static_Rotation_get_head(size_t index)
{  // ��ȡ����
	int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_GetHeadParam* param = (PT_GetHeadParam*)malloc(THREAD_NUM * sizeof(PT_GetHeadParam));  // �洢�̲߳���
	void** ret = new void* [THREAD_NUM];
	sem_init(&sem_gethead, 0, 0);
	if (tasknum < THREAD_NUM)  // ���һ���̵߳����������������߳������࣬��ô��ֱ�Ӳ������߳�
	{
		for (int i = 0; i < eks.size(); i++)  // R[lp(E[i])
		{
			if (index == eks[i].find_first())
				return eks[i];
		}
		return dynamic_bitset<>(0);
	}
	else {
		for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		{
			param[t_id].t_id = t_id;
			param[t_id].index = index;
		}
		//�����߳�
		for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		{
			pthread_create(&handles[t_id], nullptr, PT_Rotation_gethead_Inthread, (void*)&param[t_id]);
		}
		for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		{
			pthread_join(handles[t_id], &ret[t_id]);
		}
		//���߳�˯�ߣ��ȴ����еĹ����߳���ɴ�������
		for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		{
			sem_wait(&sem_gethead);
		}
		//�����ѣ�������ֵ
		for (int t_id = 0; t_id < THREAD_NUM; t_id++)
		{
			if (ret[t_id]) return *(dynamic_bitset<>*)(ret[t_id]);
		}
		return dynamic_bitset<>(0);
	}

}
//�㷨����
void GrobnerGE()
{
	dynamic_bitset<> temp = dynamic_bitset<>(n);
	for (int i = 0; i < et_num; i++)
	{
		while (ets[i].any())
		{
			temp = get_head(ets[i].find_first());
			if (!(temp.size() == 0))
				ets[i] ^= temp;
			else
			{
				eks.push_back(ets[i]);
				break;
			}
		}
	}
}
void GrobnerGE_LP()
{
	dynamic_bitset<> temp = dynamic_bitset<>(n);
	for (int i = 0; i < et_num; i++)
	{
		while (ets[i].any())
		{
			temp = get_head_lp(lp_ets[i]);  // ��Ӧ��Ԫ��
			if (temp.size())
			{
				ets[i] ^= temp;
				lp_ets[i] = ets[i].find_first();
			}
			else
			{
				eks.push_back(ets[i]);
				lp_eks.push_back(lp_ets[i]);
				lp_ets[i] = -1;
				break;
			}
		}
	}
}
void GrobnerGE_LP_static()  // ��̬λ���µĸ�˹��ȥ
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets_c[i].any())
		{
			lp = lp_ets_c[i];
			if (eks_c[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				ets_c[i] ^= eks_c[lp];
				lp_ets_c[i] = find_first_bitset(ets_c[i]);
			}
			else
			{
				eks_c[lp] = ets_c[i];
				lp_ets_c[i] = -1;
				break;
			}
		}
	}
}
void GrobnerGE_LP_static_headopt()  // ��̬λ���µĸ�˹��ȥ��ѡȡ�����һ���Ż�
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets_c[i].any())
		{
			lp = lp_ets_c[i];
			if (eks_c[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				ets_c[i] ^= eks_c[lp];
				for (int j = lp; j < column_num_c; j++)
					if (ets_c[i].test(j)) {
						lp_ets_c[i] = j;
						break;
					}
			}
			else
			{
				eks_c[lp] = ets_c[i];
				lp_ets_c[i] = -1;
				break;
			}
		}
	}
}
void GrobnerGE_LP_static_headopt_AVX()  // ��̬λ���µĸ�˹��ȥ��ѡȡ�����һ���Ż���AVX�汾
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets_c[i].any())
		{
			lp = lp_ets_c[i];
			if (eks_c[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				xor_static_AVX(ets_c[i], eks_c[lp]);
				//ets_c[i] ^= eks_c[lp];
				for (int j = lp; j < column_num_c; j++)
					if (ets_c[i].test(j)) {
						lp_ets_c[i] = j;
						break;
					}
			}
			else
			{
				eks_c[lp] = ets_c[i];
				lp_ets_c[i] = -1;
				break;
			}
		}
	}
}
void GrobnerGE_Rearrange()
{
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks_c[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{  // ����̸߳���������Ԫ�Ӳ���
			for (size_t j = 0; j < et_num_c; j++)
			{
				if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
				{
					eks_c[i] = ets_c[j];
					lp_ets_c[j] = -1;
					break;
				}
			}
		}
		for (size_t j = 0; j < et_num_c; j++)  // ѭ�����ֲ��л�
		{
			if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets_c[j] ^= eks_c[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets_c[j].test(k))
					{
						lp_ets_c[j] = k;
						break;
					}
					if (k == column_num_c - 1) lp_ets_c[j] = -1;
				}
			}
		}
	}
}
void getheadOpt_GrobnerGE()
{
	dynamic_bitset<> temp = dynamic_bitset<>(n);
	for (int i = 0; i < et_num; i++)
	{
		while (ets[i].any())
		{
			temp = PT_Static_Rotation_get_head(ets[i].find_first());  // R[lp(E[i])]������Ԫ�����뱻��Ԫ��i��������ͬ����
			if (!(temp.size() == 0))
				ets[i] ^= temp;
			else
			{
				eks.push_back(ets[i]);
				break;
			}
		}
	}
}
void test_GE()
{
	bool find_flag = false;  // true˵����������
	dynamic_bitset<>* temp = nullptr;
	for (int i = 0; i < et_num; i++)
	{
		while (ets[i].any())
		{
			find_flag = false;
			size_t index = ets[i].find_first();
			for (int j = 0; j < eks.size(); j++)  // R[lp(E[i])
			{
				if (index == eks[j].find_first())
				{
					ets[i] ^= eks[j];
					find_flag = true;
					break;
				}
			}
			if (!find_flag)
			{
				eks.push_back(ets[i]);
				break;
			}
		}
	}
}
void GrobnerGE_PT_SET()
{
	int lp = -1;
	while (!ek_set.empty())  // �������û�����������Ԫ�������
	{
		for (size_t i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��д���
		{
			if (eks_c[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
			{
				if (ek_set.find(i) != ek_set.end())
				{
					for (size_t j = 0; j < et_num_c; j++)
					{
						if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
						{
							ets_c[j] ^= eks_c[i];
							for (int k = i; k < column_num_c; k++)
								if (ets_c[j].test(k)) {
									lp_ets_c[j] = k;
									break;
								}
						}
					}
					ek_set.erase(i);  //��Ԫ�ӱ����궪��
				}
			}
			else
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (lp_ets_c[j] == i && ek_set.find(i) == ek_set.end())  // ��ֹ��Ԫ���ظ�
					{
						eks_c[i] = ets_c[j];
						lp_ets_c[j] = -1;
						ek_set.insert(i);  // ����
						break;
					}
				}
			}
		}
	}
}
void GrobnerGE_OMP()
{
	int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
	dynamic_bitset<> temp = dynamic_bitset<>(n);
	for (int i = 0; i < et_num; i++)
	{
		while (ets[i].any())
		{
			count1++;
			temp = get_head(ets[i].find_first());
			if (!(temp.size() == 0))
			{
				ets[i] ^= temp;
				count2++;
			}
			else
			{
				count3++;
				eks.push_back(ets[i]);
				break;
			}
		}
	}
	cout << count1 << " " << count2 << " " << count3 << endl;
}
void GrobnerGE_Rearrange_AVX() {
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{  // ����̸߳���������Ԫ�Ӳ���
			for (size_t j = 0; j < et_num_c; j++)
			{
				if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
				{
					eks[i] = ets[j];
					ets[j].head = -1;
					break;
				}
			}
		}
		for (size_t j = 0; j < et_num_c; j++)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				//ets[j] ^= eks[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
}
//���ݶ�ȡ����
void readData_bitset()
{
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // ��Ԫ��
	ifstream inElimTar(dir + "elimtar.txt");  // ����Ԫ��
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // �������ݶ���
	while (true)  // ��ȡ��Ԫ��
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			eks[p_ek].set(inek_loc);
		}
		if (inek.empty())
			break;
		p_ek++;
	}
	while (true)  // ��ȡ����Ԫ��
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			ets[p_et].set(inet_loc);
		}
		if (inet.empty())
			break;
		p_et++;
	}
}
void readData_reverse_bitset()
{  // �����������
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // ��Ԫ��
	ifstream inElimTar(dir + "elimtar.txt");  // ����Ԫ��
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // �������ݶ���
	while (true)  // ��ȡ��Ԫ��
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			eks[p_ek].set(n - inek_loc - 1);
		}
		if (inek.empty())
			break;
		p_ek++;
	}
	while (true)  // ��ȡ����Ԫ��
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			ets[p_et].set(n - inet_loc - 1);
		}
		if (inet.empty())
			break;
		p_et++;
	}
}
// ����������ݣ����뾲̬λ��
void readData_reverse_bitset_c_pt()
{
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // ��Ԫ��
	ifstream inElimTar(dir + "elimtar.txt");  // ����Ԫ��
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // �������ݶ���
	int lp = -1;
	while (true)  // ��ȡ��Ԫ��
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			if (lp == -1)
			{
				lp = n - inek_loc - 1;
				ek_set.insert(lp);
			}
			eks_c[lp].set(n - inek_loc - 1);
		}
		lp = -1, p_ek++;
		if (inek.empty())
			break;
	}
	while (true)  // ��ȡ����Ԫ��
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			if (lp == -1)
			{
				lp = n - inet_loc - 1;
				lp_ets_c[p_et] = lp;
			}
			ets_c[p_et].set(n - inet_loc - 1);
		}
		if (inet.empty())
			break;
		lp = -1;
		p_et++;
	}
	inElimKey.close();
	inElimTar.close();
	cout << "init_complete" << endl;
}
void init()
{
	eks.clear(); eks.resize(0);
	ets = new dynamic_bitset<>[et_num];
	lp_eks.clear(); lp_eks.resize(0);
	lp_ets = new int[et_num];
	for (int i = 0; i < ek_num; i++)
		eks.push_back(dynamic_bitset<>(n));
	for (int i = 0; i < et_num; i++)
		ets[i] = dynamic_bitset<>(n);
	readData_reverse_bitset();  // �����ʼ����Ԫ�Ӻͱ���Ԫ������
	lp_ets = new int[et_num];
	for (int i = 0; i < et_num; i++)
		lp_ets[i] = ets[i].find_first();//���㱻��Ԫ������
	for (int i = 0; i < ek_num; i++)
		lp_eks.push_back(eks[i].find_first());
}
void init_c_pt()
{
	for (int i = 0; i < column_num_c; i++)  // ��ʼ����������һ���㷨�еĲ�������
	{
		eks_c[i] = *(new bitset<column_num_c>);
		lp_eks_c[i] = -1;
	}
	readData_reverse_bitset_c_pt();  // �����ʼ����Ԫ�Ӻͱ���Ԫ������
	cout << "init_complete" << endl;
}
int main()
{
	ifstream inParam(dir + "param.txt");
	inParam >> n;  // ��������С
	inParam >> ek_num;  // ���������Ԫ�Ӹ���
	inParam >> et_num;  // ���뱻��Ԫ������
	//��̬�㷨����Ҫ�ֶ����ĳ����е�ȫ�ֳ���
	cout << "�����СΪ" << n << "����Ԫ�Ӹ���Ϊ" << ek_num << "������Ԫ������Ϊ" << et_num << endl;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

	init();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE();  // ��̬λ���洢�ľ���������˹��ȥ
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "db_GrobnerGE: " << (tail - head) * 1000.0 / freq << "ms" << endl;

	init();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_LP();  // ��̬λ���洢�ľ���������˹��ȥ��ѡȡ������Ż�
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_LP: " << (tail - head) * 1000.0 / freq << "ms" << endl;



	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_LP_static();  // ��̬λ���洢�ľ���������˹��ȥ��ѡȡ������Ż�
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_LP_static: " << (tail - head) * 1000.0 / freq << "ms" << endl;



	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_LP_static_headopt();  // ��̬λ���洢�ľ���������˹��ȥ��ѡȡ����ָ߶��Ż�
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_LP_static_headopt: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_LP_static_headopt_AVX();  // ��̬λ���洢�ľ���������˹��ȥ��ѡȡ����ָ߶��Ż�
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_LP_static_headopt_AVX: " << (tail - head) * 1000.0 / freq << "ms" << endl;



	init_c_pt();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_PT_SET();  // ��̬λ���洢�ľ���������˹��ȥ��ѡȡ����ָ߶��Ż�
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_PT_SET: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_Rearrange();  // GGELSH���߳̿黮�ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_Rearrange: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static();  // GGELSH���߳�ѭ�����ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_Block();  // GGELSH���߳̿黮�ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_Block: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE();  // ��̬λ���洢�ľ���������˹��ȥ
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE: " << (tail - head) * 1000.0 / freq<< " ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_Rearrange();  // �����ݽṹ�������㷨
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_Rearrange: " << (tail - head) * 1000.0 / freq<< " ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_Rearrange_AVX();  // �����ݽṹ�������㷨+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_Rearrange_AVX: " << (tail - head) * 1000.0 / freq<< " ms" << endl;



	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static();  // GGELSH���߳�ѭ�����ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static: " << (tail - head) * 1000.0 / freq<< " ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_Block();  // GGELSH���߳̿黮�ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_Block: " << (tail - head) * 1000.0 / freq<< " ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_AVX();  // GGELSH���߳̿黮�ְ汾+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_AVX: " << (tail - head) * 1000.0 / freq<< " ms" << endl;


	init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_Block_AVX();  // GGELSH���߳̿黮�ְ汾+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_Block_AVX: " << (tail - head) * 1000.0 / freq<< " ms" << endl;
}