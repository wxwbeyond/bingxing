#include"mybitset.h"

using namespace std;


int main(int argc, char* argv[])
{
    //静态算法则需要手动更改程序中的全局常量
    cout << "矩阵大小为" << column_num_c << "，消元子个数为" << ek_num_c << "，被消元行行数为" << et_num_c << endl;
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
    for (int i = 0; i < column_num_c; i++) // 取每个消元子，对被消元行进行操作，便于并行化
    {
        if (!eks[i].test(i)) // 消元子被逆序初始化时满足“行号” = “首项”的条件
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) // 零号进程负责升格
            {
#pragma omp single
                for (size_t j = 0; j < et_num_c; j++)
                {
                    if (i == ets[j].head)  // 说明存在对应被消元行
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
            if (rank && upshift)  // 此刻其他进程在接收数据
            {
                MPI_Recv(&(eks[i].content[0]), bit_size, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(ets[i].head), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            upshift = 0;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank)
        {  // 0号进程以外的进程做消去
            cout << "rank = " << rank << " come2" << endl;
#pragma omp for
            for (int j = rank - 1; j < et_num_c; j += PROGRESS_NUM - 1) // 循环划分并行化处理被消元行
            {
                if (i == ets[j].head)  // 说明存在对应被消元行
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
                    //将计算后的结果广播回0号进程
                    MPI_Send(&(ets[j].content[0]), bit_size, MPI_INT, 0, j, MPI_COMM_WORLD);
                    MPI_Send(&(ets[j].head), 1, MPI_INT, 0, j, MPI_COMM_WORLD);
                }
            }
        }
        else 
        {  // 结果存在0号进程
#pragma omp for
            for (int j = 0; j < et_num_c; j++) // 处理数据
            {
                if (i == ets[j].head)  // 说明存在对应被消元行并且子进程执行了操作
                {
                    //获取计算后的结果
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

