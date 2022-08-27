from mimetypes import init
import pandas as pd
import numpy as np
import sys
import mykmeanssp



def calc_distance(vec1, vec2):#calculate distance between two vectors 
    return sum((vec1 - vec2)**2)
          
def kmeans_pp(k,data_arr,n,d): #init the first k centroids
    np.random.seed(0)
    first_centroid = np.random.choice(n)

    centroids_index = []
    p_list =  np.zeros((n), dtype=np.float64) 
    d_list =  np.full(n, np.inf) 

    centroids_index.append(first_centroid)
    d_sum = 0
    j = 1
    while j < k:
        new_centroid = data_arr.loc[centroids_index[j-1]]
        for i in range(n):
            curr_distance = calc_distance(data_arr.loc[i], new_centroid)
            if curr_distance < d_list[i]:
                d_list[i] = curr_distance
        d_sum = d_list.sum()

        for i in range(n):
            p_list[i] = d_list[i] / d_sum
        next_centroid = np.random.choice(n, size=None, p=p_list)
        centroids_index.append(next_centroid)
        j += 1
    
    return centroids_index, data_arr.loc[centroids_index]

def recieve_input():#use the system atguments as variables 
    arguments_size = len(sys.argv)
    k_float = float(sys.argv[1])
    k = int(k_float)   
    if k_float != k:
        print("Invalid Input")
        exit()
    if (arguments_size == 6):
        epsilon = float(sys.argv[3])
        max_iter = int(sys.argv[2])
        input_file_1 = sys.argv[4]
        input_file_2 = sys.argv[5]
    elif (arguments_size == 5):
        epsilon = float(sys.argv[2])
        max_iter = 300
        input_file_1 = sys.argv[3]
        input_file_2 = sys.argv[4]
    else:
        print("Invalid Input")
        exit()
    
    return k, max_iter,epsilon, input_file_1, input_file_2

def validate_input(k,max_iter,n,epsilon): #validate the input 
    if(k < 1 or k >= n or max_iter < 0):
        print("Invalid Input")
        exit()
    return

def check_len(coordinate):# using for printing format
            str_coor = str(coordinate)
            str_coor_array = str_coor.split(".")
            curr_len = len(str_coor_array[1])
            if curr_len < 4:
                str_coor = str_coor + '0'    
            return(str_coor)

def print_final_centroids(k_centriods):#print the final k centroids acording to the format
    for centroid in k_centriods:
        print(",".join([check_len(round(coordinate,4)) for coordinate in centroid]))

if __name__ == "__main__":
    k, max_iter,epsilon, input_file_1, input_file_2 = recieve_input() 
    data_frame1 = pd.read_csv(input_file_1 , header=None, index_col=0)
    data_frame2 = pd.read_csv(input_file_2 , header=None, index_col=0)
    merge_data_frame = pd.merge(data_frame1, data_frame2, left_index=True, right_index=True)
    n = len(merge_data_frame) # n = num of rows
    d = len(merge_data_frame.columns)# d=  num of columns
    validate_input(k,max_iter,n,epsilon)
    first_init, first_k = kmeans_pp(k,merge_data_frame, n,d)
    vectors = merge_data_frame.values.tolist() # DataFrame -> [[],[],[]...]
    k_init_centroids = first_k.values.tolist() 
    k_centroids = mykmeanssp.fit(n, d, k, max_iter, epsilon, vectors, k_init_centroids)
    print(",".join([str(i) for i in first_init]))
    print_final_centroids(k_centroids)
    


    

    

    
    


    
    





