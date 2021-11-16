from operator import le
from re import T
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import partition, size
from numpy.lib.polynomial import RankWarning
from numpy.lib.twodim_base import mask_indices

 
PRINT_MAXWELL = False
PRINT_SRED_R_2 = False
PRINT_ENERGY = False
PRINT_RADIAL = False
PRINT_AKFC = False

R_size_sosyd = 5
kolvo_mol = 64
delta_time = 0.001



kakoi_print = 2
if(kakoi_print != 0):
    if(kakoi_print == 1):
        PRINT_MAXWELL = True
    elif(kakoi_print == 2):
        PRINT_SRED_R_2 = True
    elif(kakoi_print == 3):
        PRINT_ENERGY = True
    elif(kakoi_print == 4):
        PRINT_RADIAL = True
    elif(kakoi_print == 5):
        PRINT_AKFC = True

#data = np.loadtxt('C:/Users/Компьютер/Desktop/total_energy.txt', delimiter='\t', dtype=np.double)



#data_Ek_Ep = np.loadtxt('C:/code_bocks/программы/gorizonti/Data/energy.txt', delimiter='\t', dtype=np.double)
#Maxwell_source = 'C:/code_bocks/программы/gorizonti/Data/energy.txt'

def get_maxwell_source(n):
    return 'C:/code_bocks/программы/gorizonti/Data/Maxwell/Maxwell_' + str(int(n)) + '.txt'


def get_radial_source(n):
    return 'C:/code_bocks/программы/gorizonti/Data/Radial_function/Radial_' + str(int(n)) + '.txt'






def split(itog_V, probnie_V, V_max, probnie_V_max, lenght):
    delta_V = V_max / lenght
    probnie_delta_V = probnie_V_max / lenght 

    i = 0 #нумеровка в итоговом массиве
    j = 0
    do_this = True
    while(do_this):
        if((j+1)*probnie_delta_V <= (i+1)*delta_V):
            itog_V[i] += probnie_V[j]
            j += 1

        else:
            l_border = j*probnie_delta_V
            while((j+1)*probnie_delta_V > (i+1)*delta_V):
                itog_V[i] += ((i+1)*delta_V - l_border)/probnie_delta_V * probnie_V[j]
                l_border = (i+1)*delta_V
                i += 1

                if(i == lenght):
                    do_this = False
                    break

            if(i < lenght):
                itog_V[i] += ((j+1)*probnie_delta_V - l_border)/probnie_delta_V * probnie_V[j]
                j += 1


        if(j == lenght):
            do_this = False



#---------------------[НАЧАЛЬНЫЕ ПАРАМЕТРЫ, ДЛЯ МАКСВЕЛЛА]-----------------------------
start = 700 #включительно
end = 850   #включительно

counter = end - start + 1 #количество наборов


source_first = open(get_maxwell_source(1)) #на самом деле можно любой открыть, так как max_R и lenght везде одинаковые

lenght = -1 #сколько в одном наборе
for line in source_first:
    lenght += 1


data_V_max = [0] * counter
V_max = 0
data_V_for_hist = [[0] * lenght for i in range(counter)]

#-------------------------------------------------------------------------------------
for j in range(start, end + 1):
    l = 0
    source = open(get_maxwell_source(j))
    for line in source:
        if(l == 0):
            data_V_max[j - start] = np.double(line) 
            V_max += data_V_max[j - start]
        else:
            data_V_for_hist[j - start][l-1] = int(line)

        l += 1


V_max = V_max / counter
#----------------------------[ВЫВОД МАКСВЕЛЛА]----------------------------------------
data_total_hist = [0] * lenght

for j in range(counter):
    split(data_total_hist, data_V_for_hist[j], V_max, data_V_max[j], lenght)


if(PRINT_MAXWELL):
    for i in range(0, lenght):
        plt.plot(i*V_max/lenght, data_total_hist[i]/counter, 'go')
#-------------------------------------------------------------------------------------







#------------------------------[СРЕДНЕЕ R^2]------------------------------------------
if(PRINT_SRED_R_2):
    data_sred_r_2 = np.loadtxt('C:/code_bocks/программы/gorizonti/Data/Sred_r_2.txt', delimiter='\t', dtype=np.double)

    for i in range(0, len(data_sred_r_2)):
        plt.plot(i * delta_time, data_sred_r_2[i], 'ro')
        plt.plot(i * delta_time, 6.4 * i * delta_time, 'go')
#-------------------------------------------------------------------------------------





#-------------------------------[ЭНЕРГИЯ]---------------------------------------------
if(PRINT_ENERGY):
    data_energy = np.loadtxt('C:/code_bocks/программы/gorizonti/Data/energy.txt', delimiter='\t', dtype=np.double)

    for i in range(0, len(data_energy), 2):
        plt.plot(i * delta_time, data_energy[i] + data_energy[i+1], 'ro')

#-------------------------------------------------------------------------------------



#------------------------------[РАДИАЛЬНАЯ ФУНКЦИЯ РАСПРЕДЕЛЕНИЯ]---------------------
start = 200 #включительно
end = 500   #включительно

if(PRINT_RADIAL):
    counter = end - start + 1 #количество наборов

    source_first = open(get_radial_source(1)) #на самом деле можно любой открыть, так как max_R и lenght везде одинаковые

    max_R = 0
    lenght = -1 #сколько в одном наборе

    for line in source_first:
        if(lenght == -1):
            max_R = np.double(line) 


        lenght += 1

    data_n_for_radial = [0] * (lenght+1)


    for j in range(start, end + 1):
        l = 0
        source = open(get_radial_source(j))
        for line in source:
            if(l != 0):
                data_n_for_radial[l-1] += int(line)

            l += 1

    
    data_n_it = [0] * lenght
    data_n_it[0] = data_n_for_radial[0] + data_n_for_radial[1]
    for i in range(1,lenght):
            data_n_it[i] = data_n_it[i-1] + data_n_for_radial[i+1]
            #print(data_n_it[i-1], data_n_it[i])


    for i in range(0, lenght-1):
        plt.plot((i+1)/lenght * max_R, data_n_it[i] / counter / (4/3 * 3.14 * kolvo_mol * ((i+1)/lenght * max_R / R_size_sosyd)**3) / kolvo_mol, 'ro') 
#-------------------------------------------------------------------------------------



#------------------------------[АКФС]-------------------------------------------------
if(PRINT_AKFC):
    data_akfc= np.loadtxt('C:/code_bocks/программы/gorizonti/Data/AKFC.txt', delimiter='\t', dtype=np.double)

    integ = 0

    for i in range(0, len(data_akfc)):
        plt.plot(i * delta_time, data_akfc[i], 'ro')
        integ += data_akfc[i] * delta_time
        
    print(integ)

    

#-------------------------------------------------------------------------------------



plt.show()
