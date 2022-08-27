from mimetypes import init
import pandas as pd
import numpy as np
import sys
import mykmeanssp



def initial_Points(filename):
     vectors = []
     file = open(filename, 'r')
     for row in file:
         filter_Row = row.split(',')
         new_Row=[]
         for p in filter_Row:
             new_Row.append(float(p))
         vectors.append(new_Row)
     file.close()
     return vectors

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
    if arguments_size ==4:
        k_float = float(sys.argv[1])
        k = int(k_float)   
        if k_float != k:
            print("Invalid Input")
            exit()   
        flow = sys.argv[2]
        input = sys.argv[3]      

    elif arguments_size ==3:
        k=0
        flow = sys.argv[1]
        input = sys.argv[2]

    else :
        print("Invalid Input!")
        exit()
    
    return k, flow,input

def validate_input(k): #validate the input 
    if(k < 1 or k >= n):
        print("Invalid Input")
        exit()
    return

if __name__ == "__main__":
    try:
        vectors = initial_Points(input)
    except:
        print("An Error Has Occured")
        exit()
    n = len(vectors) # n = num of rows
    d = len(vectors[0])# d=  num of columns
    k,flow, input = recieve_input() 
    validate_input(k,n)
    max_iter=300

    try:
        if flow=="spk":
            newDataPoints = spkmeans.getnew_datapoints(n,d,k,vectors)
            data_Points_df = pd.DataFrame(newDataPoints)
            k=len(newDataPoints[0])
            centroids_indexes,centroids_points  = kmeans_pp(k, data_Points_df,n,d)
            centroids_list =centroids_points.values.tolist()
            for i in range(k - 1):
                print(centroids_indexes[i], end=',')
            print(centroids_indexes[k - 1])
            k_centroids = spkmeans.kmeans_pp(n, k, k, max_iter, newDataPoints, centroids_list)
        elif flow=="wam":
            spkmeans.wam(n,d,vectors)
        elif flow=="ddg":
            spkmeans.ddg(n,d,vectors)
        elif flow == "lnorm":
            spkmeans.lnorm(n,d,vectors)
        elif flow=="jacobi":
            spkmeans.jacobi(n,vectors)
    except:
        print("An Error Has Occurred")
        exit()   
