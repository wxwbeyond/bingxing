#pragma warning(disable : 4996)
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include"mybitset.h"
using namespace std;

const int BLOCK_SIZE = 1024;
int** eks_bits = new int* [column_num_c];
int** ets_bits = new int* [et_num_c];
//CUDA�㷨�˺���
__global__ void upshift_kernel(int i, int** gpu_eks, int** gpu_ets, int* gpu_lp_ets, int et_num_c) 
{
	//���߳���
	if (blockIdx.x == 0 && threadIdx.x == 0)
		for (int j = 0; j < et_num_c; j++)
		{
			if (i == gpu_lp_ets[j])  // ˵�����ڶ�Ӧ����Ԫ��
			{
				gpu_eks[i] = gpu_ets[j];
				gpu_lp_ets[j] = -1;
				break;
			}
		}
}
__global__ void elim_kernel(int i, int** gpu_ek, int** gpu_et, int* gpu_lp_ets, int et_num_c, int column_num_c, int bit_size)
{
	int tx = blockDim.x * blockIdx.x + threadIdx.x;
	int row = blockIdx.x;//ÿ���鸺��һ������Ԫ��
	bool find = false;
	for (int j = row; j < et_num_c; j += blockDim.x)  // ѭ�����ֲ��л�
	{
		if (i == gpu_lp_ets[j])  // ˵�����ڶ�Ӧ����Ԫ��
		{
			for (int k = 0; k < bit_size; k++)
				gpu_et[j][k] ^= gpu_ek[i][k];
			__syncthreads();
			if (threadIdx.x == 0)
			{
				for (int k = 0; k < bit_size; k++)
				{
					for (int l = 0; l < 32; l++)
					{
						if ((gpu_et[j][k] & (1 << l)))
						{
							gpu_lp_ets[j] = k * 32 + l;
							find = true;
							break;
						}
					}
					if (find)
					{
						find ^= find;
						break;
					}
					else
						gpu_lp_ets[j] = -1;
				}
			}
		}
	}
}
//CUDA��ȥ�㷨
void CUDA_GE(bitset<column_num_c>* ek, bitset<column_num_c>* et, int* lp_ets)
{
	for (int i = 0; i < column_num_c; i++)
	{
		eks_bits[i] = new int[bit_size];
		for (int j = 0; j < bit_size; j++)
			eks_bits[i][j] = 0;
	}
	for (int i = 0; i < et_num_c; i++)
	{
		ets_bits[i] = new int[bit_size];
		for (int j = 0; j < bit_size; j++)
			ets_bits[i][j] = 0;
	}
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
				lp = column_num_c - inek_loc - 1;
				lp_eks_c[lp] = lp;
			}
			// ���set����
			eks_bits[lp][(column_num_c - inek_loc - 1) / 32] |= (1 << ((column_num_c - inek_loc - 1) % 32));
		}
		lp = -1;
		p_ek++;
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
				lp = column_num_c - inet_loc - 1;
				lp_ets_c[p_et] = lp;
			}
			ets_bits[p_et][(column_num_c - inek_loc - 1) / 32] |= (1 << ((column_num_c - inek_loc - 1) % 32));
		}
		lp = -1;
		p_et++;
		if (inet.empty())
			break;
	}
	inElimKey.close();
	inElimTar.close();
	cudaError_t ret;//���ڴ����飬�� CUDA �ӿڵ��óɹ��᷵�� cudaSucess
	int** gpu_ets;
	int** gpu_eks;
	int** gpu_ets_bits = new int* [et_num_c];  // ��ŵ���ָ���׵�ַ
	int** gpu_eks_bits = new int* [column_num_c];  // ��ŵ���ָ���׵�ַ
	int* gpu_lp_ets;
	int gpu_et_size = et_num_c * sizeof(int*);
	int gpu_ek_size = column_num_c * sizeof(int*);
	int bitset_size = sizeof(int) * bit_size;
	int gpu_lp_ets_size = et_num_c * sizeof(int);
	for (int i = 0; i < et_num; i++)
	{
		int* bits;  // ��ʱ����Դ��ַ
		if (cudaMalloc(&bits, gpu_et_size) != cudaSuccess)
			printf("cudaMalloc gpudata failed!\n");
		if (cudaMemcpy(bits, ets_bits[i], bitset_size, cudaMemcpyHostToDevice) != cudaSuccess)
			printf("cudaMemcpyHostToDevice failed!\n");
		gpu_ets_bits[i] = bits;
	}
	if (cudaMalloc(&gpu_ets, gpu_et_size) != cudaSuccess)
		printf("cudaMalloc gpudata failed!\n");
	if (cudaMemcpy(gpu_ets, gpu_ets_bits, gpu_et_size, cudaMemcpyHostToDevice) != cudaSuccess)
		printf("cudaMemcpyHostToDevice failed!\n");
	for (int i = 0; i < column_num_c; i++)
	{
		int* bits;  // ��ʱ����Դ��ַ
		if (cudaMalloc(&bits, gpu_ek_size) != cudaSuccess)
			printf("cudaMalloc gpudata failed!\n");
		if (cudaMemcpy(bits, eks_bits[i], bitset_size, cudaMemcpyHostToDevice) != cudaSuccess)
			printf("cudaMemcpyHostToDevice failed!\n");
		gpu_eks_bits[i] = bits;
	}
	if (cudaMalloc(&gpu_eks, gpu_ek_size) != cudaSuccess)
		printf("cudaMalloc gpudata failed!\n");
	if (cudaMemcpy(gpu_eks, gpu_eks_bits, gpu_ek_size, cudaMemcpyHostToDevice) != cudaSuccess)
		printf("cudaMemcpyHostToDevice failed!\n");
	//�����Դ�ռ䲢�ҽ��д�����
	if (cudaMalloc(&gpu_lp_ets, gpu_lp_ets_size) != cudaSuccess)
		printf("cudaMalloc gpudata failed!\n");
	//�����ݴ����� GPU �˲����д�����
	if (cudaMemcpy(gpu_lp_ets, lp_ets_c, gpu_lp_ets_size, cudaMemcpyHostToDevice) != cudaSuccess)
		printf("cudaMemcpyHostToDevice failed!\n");
	dim3 dimBlock(BLOCK_SIZE, 1), dimGrid(1, 1); //�߳̿顢�߳�����
	cudaEvent_t start, stop;  //��ʱ��
	float elapsedTime = 0.0;
	cudaEventCreate(&start), cudaEventCreate(&stop);
	cudaEventRecord(start, 0);  //��ʼ��ʱ
	cudaError_t exec;
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!(eks_bits[i][i / 32] & (1 << (i % 32))))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
			upshift_kernel << <1, 1 >> > (i, gpu_eks, gpu_ets, gpu_lp_ets, et_num_c);
		cudaDeviceSynchronize();//CPU �� GPU ֮���ͬ������
		exec = cudaGetLastError();
		if (exec != cudaSuccess)
			printf("upshift_kernel failed, %s\n", cudaGetErrorString(exec));
		elim_kernel << <dimGrid, dimBlock >> > (i, gpu_eks, gpu_ets, gpu_lp_ets, et_num_c, column_num_c, bit_size);//������ȥ����ĺ˺���
		cudaDeviceSynchronize();//CPU �� GPU ֮���ͬ������
		exec = cudaGetLastError();
		if (exec != cudaSuccess)
			printf("elim_kernel failed, %s\n", cudaGetErrorString(exec));
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);//ֹͣ��ʱ
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("CUDA_GE:%f ms\n", elapsedTime);
	cudaError_t cudaStatus2 = cudaGetLastError();
	if (cudaGetLastError() != cudaSuccess)
		fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus2));
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
}
int main()
{
	cout << "�����СΪ" << column_num_c << "����Ԫ�Ӹ���Ϊ" << ek_num_c << "������Ԫ������Ϊ" << et_num_c << endl;
	CUDA_GE(eks_c, ets_c, lp_ets_c);
}