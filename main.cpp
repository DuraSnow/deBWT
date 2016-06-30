//
//  main.cpp
//  deBWT
//
//  Created by Abelard on 16/6/20.
//  Copyright © 2016年 Abelard. All rights reserved.
//

#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <bitset>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
using namespace std;
size_t KMER=20;
size_t KMER_N=4864385;
uint64_t *io_K;
size_t io_index = 0,bc_index;
string branch_c = "";
uint64_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    3};

bool cmp(const uint64_t a, const uint64_t b)
{
   
    return  (a<<2) < (b<<2);
}
bool bwt_cmp(const uint64_t a, const uint64_t b)
{
    uint64_t ta=(a>>3), tb=(b>>3);
    while (ta<bc_index && tb<bc_index && branch_c[ta] == branch_c[tb])
    {
        ++ta;
        ++tb;
    }
    if (ta==bc_index) return true;
    if (tb==bc_index) return false;
    return (branch_c[ta] < branch_c[tb]);
}

size_t search_k(uint64_t n)
{
    //binary search
    size_t l=0, r=io_index-1, m;
    while (l < r)
    {
        m = (l+r)/2;
        if (n <= io_K[m]) r = m;
        else l = m + 1;
    }
    if (io_K[r] == n) return r;
    return -1;
}
int main()
{
    ifstream f;
    f.open("/Users/abelard/Desktop/E.coli.fa");
    string buffer,dnaref;
    getline(f,buffer);
    while(getline(f,buffer)) dnaref +=buffer;
   
    f.close();
    size_t ref_l = dnaref.length();
    FILE *f_kmer;
    f_kmer = fopen("/Users/abelard/Desktop/mer_counts_dumps.fa","rb");
    char curkmer[100];
    size_t klen= 0,j;
    
    uint64_t *all_K,kmerindex;
    uint64_t tmp = 0, cur;
    all_K = new uint64_t[ref_l];
    
    while (fscanf(f_kmer, ">%lld\n%s\n", &kmerindex, curkmer) ==2)
    {
        tmp=0;
        j = 0;
        for (size_t j=0; j<KMER+2; ++j)tmp = (tmp << 2) | get_c[curkmer[j]];
        cur = tmp << 20;
        all_K[klen]=cur|(kmerindex);
        
        klen++;
    }
    KMER_N = klen;
    sort(all_K,all_K+klen,cmp);
    
    uint64_t mask_k=-1, mask_out=-1, mask_in=-1, mask_n=-1;
    mask_k = mask_k << 2 >> 2 >> (64 - 2*(KMER+2) + 2) << (64 - (KMER+2)*2 + 2);
    mask_in = mask_in >> 62 << 62;
    mask_out = mask_out >> (64 - 2*(KMER+2)) << (64 - 2*(KMER+2)) << (KMER+1)*2 >> (KMER+1)*2;
    mask_n = mask_n << (2*(KMER+2)) >> (2*(KMER+2));
    
    uint64_t io_cnt[100000][4],now_kmer,last_kmer;
    
    io_K= new uint64_t[1000];
    uint64_t tmp_num,Kout_index=0;
    bool Multi_in, Multi_out;
    for (size_t i=0; i<KMER; ++i)
    {
        tmp = (tmp << 2) | get_c[dnaref[i]];
    }
    uint64_t begin_kmer = tmp << (64 - 40) >> 2;

    for (size_t i=0, j, k; i<KMER_N;)
    {
        Multi_in=false; Multi_out=false;
        tmp_num=0;


        j = i+1;
        now_kmer = all_K[i]&mask_k;
        last_kmer = now_kmer;
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer)
        {
            if (now_kmer == begin_kmer)
            {
                Multi_in = true;
                if ((all_K[i]&mask_out)>>(64-(KMER+2)*2) != get_c[dnaref[KMER]]) Multi_out = true;
                ++tmp_num;
            }
            
        }
        while (j<klen && ((all_K[j]&mask_k) == (now_kmer)))
        {
            if ((all_K[j]&mask_in) != (all_K[i]&mask_in)) Multi_in = true;
            if ((all_K[j]&mask_out) != (all_K[i]&mask_out)) Multi_out = true;
            ++j;
            
        }
        
        if (Multi_in || Multi_out)
        {
      
            if (Multi_out) io_cnt[io_index][0]=1;
            if (Multi_in)
            {
                for (k = i; k<j; ++k)
                {
                    tmp_num += (all_K[k]&mask_n);
                }
                io_cnt[io_index][1]=1;
                io_cnt[io_index][2]=tmp_num;
                io_cnt[io_index][3]=Kout_index;
                Kout_index += tmp_num;
               
                
            }
            io_K[io_index] = (all_K[i]&mask_k) >> (64-(KMER+2)*2+2);
            io_index++;
            cout<<io_index<<endl;
          
        }
        i = j;
    }
    cout<<io_index<<endl;
    uint64_t *MulitinBC;
    size_t MulitinBC_size = (Kout_index);
    MulitinBC = new uint64_t[MulitinBC_size];
    for (size_t i=0; i<MulitinBC_size; ++i) MulitinBC[i]=0;
    uint64_t mask_rk = -1;
    mask_rk = mask_rk << (64-KMER*2) >> (64-KMER*2);
    bc_index=1;
    uint64_t MulitinBC_index, MulitinBC_len, tk;
    tmp = 0;
    for (size_t i=0, index; i<ref_l-1; ++i)
    {
        tmp = (tmp << 2) | get_c[dnaref[i]];
        if (i>=KMER-1)
        {
            cur = tmp&mask_rk;
            index = search_k(cur);
            if (index!=-1)
            {
                if (io_cnt[index][1])
                {
                    
                    MulitinBC_index = (io_cnt[index][3]);
                    MulitinBC_len = ((io_cnt[index][2]));
                   
                    tk = MulitinBC_index;
                 
                    while (tk<MulitinBC_size && MulitinBC[tk])
                        ++tk;
                    if (tk < (MulitinBC_index + MulitinBC_len))
                    {
                       
                        if (i>=KMER)
                            MulitinBC[tk] = (bc_index << 3)|get_c[dnaref[i-KMER]];
                        else 
                            MulitinBC[tk] = (bc_index << 3)|4;
                      
                    }
                    else
                    {
                        cout << "something wrong! attention! build MulitinBC out of range\n";
                  
                    }

                }
             
                if (io_cnt[index][0])
                {
                
                  
                    branch_c+=dnaref[i+1];
                    bc_index++;
                    
                }
            }
        }
    }

    
    uint64_t mask_code = -1;
    mask_code = mask_code << 61 >> 61;
    last_kmer = 0;

    uint8_t *BWT, tmp_code;
    string BWT_string="";
    BWT = new uint8_t[ref_l];
    size_t bwt_index = 0;

    uint64_t *last_suffix;
    last_suffix = new uint64_t[KMER+1];
    for (size_t i=KMER+1; i>=1; --i)
    {
        for (size_t j=ref_l-i; j<ref_l; ++j) tmp = (tmp << 2) | get_c[dnaref[j]];
        last_suffix[KMER+1-i] = tmp<<(64-(i)*2);
       
    }
    sort(last_suffix, last_suffix+KMER+1, cmp);

    


    for (size_t i=0, j, tmp_index=0, l_index=0; i<klen;)
    {
   
        Multi_in=false; Multi_out=false; 
        j = i+1;
        now_kmer = all_K[i]&mask_k;
        while (l_index < KMER+1 && (last_suffix[l_index]&mask_k) <= now_kmer)
        {
            BWT[bwt_index++] = ((last_suffix[l_index++]&mask_in)>>62);
            
        }
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer) 
        {
            if (now_kmer == begin_kmer)
            {
                Multi_in = true;
            }
            else
            {
                BWT[bwt_index++] = 4;
            }
          
        }
        last_kmer = now_kmer;
        while (j<klen && ((all_K[j]&mask_k) == now_kmer))
        {
            if ((all_K[j]&mask_in) != (all_K[i]&mask_in)) Multi_in = true;
            if ((all_K[j]&mask_out) !=( all_K[i]&mask_out)) Multi_out = true;
            ++j;
        }
        if (Multi_in || Multi_out)
            tmp_index++;
        if (Multi_in)
        {
            MulitinBC_index = (io_cnt[tmp_index][3]);
            MulitinBC_len = (io_cnt[tmp_index][2]);
            
            sort(MulitinBC+MulitinBC_index, MulitinBC+MulitinBC_index+MulitinBC_len, bwt_cmp);

            for (size_t k=MulitinBC_index; k<MulitinBC_index+MulitinBC_len; ++k)
            {
                BWT[bwt_index++] = MulitinBC[k]&mask_code;
            }
        }
        else
        {
            tmp_num=0;  
            for (size_t k = i; k<j; ++k)
            {
                tmp_num += (all_K[k]&mask_n);
            }
            tmp_code = ((all_K[i]&mask_in)>>62);
            while (tmp_num--)
            {
                BWT[bwt_index++] = tmp_code;
            }
        }
        i = j;
    }
    char c[6]={'A', 'C', 'G', 'T', '$', 'X'};
    for (size_t i=0; i<=ref_l; ++i)
    {
        
        BWT_string +=c[BWT[i]];
    }
    ofstream outfile("/Users/abelard/Desktop/bwt.fa",ios::out);
    outfile<<BWT_string<<endl;
    printf("Time used = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    delete [] all_K;
    delete [] io_K;
    branch_c.erase();
    delete [] BWT;
    return 0;
    
}

