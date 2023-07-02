#include"mybitset.h"

using namespace std;


int main(int argc, char* argv[])
{
    //��̬�㷨����Ҫ�ֶ����ĳ����е�ȫ�ֳ���
    cout << "�����СΪ" << column_num_c << "����Ԫ�Ӹ���Ϊ" << ek_num_c << "������Ԫ������Ϊ" << et_num_c << endl;
    init_MyB();
    int provided;
    MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
    double head, tail;
    head = MPI_Wtime();
    MPI_Request request;
    int rank = 0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int upshift = 0;
    int i, j, k;
#pragma omp parallel num_threads(THREAD_NUM), private(i, j, k)
    for (int i = 0; i < column_num_c; i++) // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
    {
        if (!eks[i].test(i)) // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) // ��Ž��̸�������
            {
#pragma omp single
                for (size_t j = 0; j < et_num_c; j++)
                {
                    if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
                    {
                        eks[i] = ets[j];
                        ets[j].head = -1;
                        upshift = 1;
                        for (int j = 1; j < PROGRESS_NUM; j++)
                        {
                            MPI_Send(&upshift, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                            MPI_Send(&(eks[i].content[0]), bit_size, MPI_INT, j, 1, MPI_COMM_WORLD);
                            MPI_Send(&(ets[i].head), 1, MPI_INT, j, 2, MPI_COMM_WORLD);
                        }
                        break;
                    }
                }
                if (!upshift)
                {
                    for (int j = 1; j < PROGRESS_NUM; j++)
                        MPI_Send(&upshift, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                }
            }
            else
                MPI_Recv(&upshift, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (rank && upshift)  // �˿����������ڽ�������
            {
                MPI_Recv(&(eks[i].content[0]), bit_size, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(ets[i].head), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            upshift = 0;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank)
        {  // 0�Ž�������Ľ�������ȥ
            cout << "rank = " << rank << " come2" << endl;
#pragma omp for
            for (int j = rank - 1; j < et_num_c; j += PROGRESS_NUM - 1) // ѭ�����ֲ��л�������Ԫ��
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
                        if (k == column_num_c - 1) 
                            ets[j].head = -1;
                    }
                    //�������Ľ���㲥��0�Ž���
                    MPI_Send(&(ets[j].content[0]), bit_size, MPI_INT, 0, j, MPI_COMM_WORLD);
                    MPI_Send(&(ets[j].head), 1, MPI_INT, 0, j, MPI_COMM_WORLD);
                }
            }
        }
        else 
        {  // �������0�Ž���
#pragma omp for
            for (int j = 0; j < et_num_c; j++) // ��������
            {
                if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ�в����ӽ���ִ���˲���
                {
                    //��ȡ�����Ľ��
                    cout << "receiving:  " << j << endl;
                    MPI_Recv(&(ets[j].content[0]), bit_size, MPI_INT, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&(ets[j].head), 1, MPI_INT, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }
    tail = MPI_Wtime();
    MPI_Finalize();
    if (!rank)
    {
        cout << "GrobnerGE_MPI: " << (tail - head) * 1000 << " ms" << std::endl;
    }
    return 0;
}

