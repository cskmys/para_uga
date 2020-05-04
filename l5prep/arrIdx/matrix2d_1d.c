
#include<stdio.h>
#include<assert.h>

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

int getIdx(int nbEle, int idx){
    idx = idx % nbEle;
    if(idx < 0){
        idx = nbEle + idx;
    }
    return idx;
}

int getEle(int *arr, int nbEle, int idx){
    idx = getIdx(nbEle, idx);
    return arr[idx];
}

int main(){
    int n = 8;
    int p = n / 4;

    assert(n % p == 0); // important simplification we are making

    int arr[n * n];
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            int idx = (i * n) + j;
            arr[idx] = idx;
            printf("%02d ", arr[idx]);
        }
        printf("\n");
    }

    int brr[n * n];
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            int idx = (i * n) + j;
            brr[idx] = 0;
            if(i == j){
                brr[idx] = 1;
            }
            printf("%02d ", arr[idx]);
        }
        printf("\n");
    }


    
    printf("-------------------------------\n");
    for(int rank = 0; rank < p; ++rank){
        int sub_array[n/p * n];
        for(int i = 0; i < n/p; ++i){
            for(int j = 0; j < n; ++j){
                int main_arr_idx = ((i + (rank * n/p)) * n) + j;
                int sub_arr_idx = (i * n) + j;
                sub_array[sub_arr_idx] = arr[main_arr_idx];
                printf("%02d ", sub_array[sub_arr_idx]);
            }
            printf("\n");
        }
        printf("\n");
        for(int step = 0; step < p; ++step){
            for (int i = 0; i < n/p; ++i){
                for (int j = 0; j < n/p; ++j){
                    int idx = getIdx( n, (((rank - step) % p) * (n/p)) + j);
                    int sub_arr_idx = (i * n) + idx;
                    printf("%02d ", sub_array[sub_arr_idx]);
                }
                printf("\n");
            }
            printf("\n");
        }
        printf("-------------------------------\n");
    }
    printf("\n");
    return 0;
}