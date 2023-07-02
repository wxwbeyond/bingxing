#pragma once
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<ctime>
#include<bitset>
#include<vector>
#include<string>
#include<sstream>
#include<unordered_set>
#include<windows.h>
#include<immintrin.h> 
#include<pthread.h>
#include<semaphore.h>
#include<mutex>
#include<condition_variable>
#include<cmath>
#include<stdio.h>
#include<time.h>
#include<algorithm>
#include<mpi.h>
#include<omp.h>
#include<boost/dynamic_bitset.hpp>
using namespace std;
#ifndef CIRCLE_H
#define CIRCLE_H

#define THREAD_NUM 8   // �߳�����
#define PROGRESS_NUM 6 // ��������

int ek_num;  // ������Ԫ�Ӹ���
int et_num;  // ���뱻��Ԫ������

int n = 0;  // �����С
const int k = 1;

const int column_num_c = 23075;
const int ek_num_c = 18748;  // ������Ԫ�Ӹ���
const int et_num_c = 14325;  // ���뱻��Ԫ������
int bit_size = column_num_c / 32 + 1;


string dir = "/data/a10/";
stringstream ss;

bitset<column_num_c> eks_c[column_num_c]; // ��Ԫ�ӣ�����һЩ���ڼ���������
bitset<column_num_c> ets_c[et_num_c];

int lp_ets_c[et_num_c];
int lp_eks_c[column_num_c];

long long head, tail, freq;

class MyBitSet
{
public:
	int head;  
	int* content;
	MyBitSet();
	MyBitSet& operator^=(const MyBitSet& b);
	MyBitSet& my_xor_AVX(const MyBitSet& b);
	int test(int index);
	void set(int index);
	bool any();
};
MyBitSet eks[column_num_c];
MyBitSet ets[et_num_c];

int find_first_bitset(const bitset<column_num_c>& b);
void readData_reverse_bitset_c();
void init_c();
void GrobnerGE_OMP();
void readData_reverse_MyB();
void init_MyB();
#endif