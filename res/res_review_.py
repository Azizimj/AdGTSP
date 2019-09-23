# import csv
#
# with open('res.csv','rb') as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')
#     headers = next(readCSV, None)
#     for row in readCSV:
#         print(row)
#         print(row[0])
#         print(row[0],row[1],row[2],)
#
# csvfile.close()

import pandas, numpy as np, re, matplotlib.pyplot as plt


def extract_plot(name, sig, ave):
    print(name+" "+str(sig)+" "+str(ave))
    chosen = np.zeros((1, dim))
    for i, row in enumerate(rd[name]):
        if str(row) != 'nan' and str(rd['Algo'][i]) == 'AdNNnew ':
            # print(row)
            tmp = re.findall('\.\d+',row)
            tmp = list(map(float, tmp))
            len_ = len(tmp)
            tmp = np.array(tmp).reshape(int(len_/dim),dim)
            chosen = np.vstack((tmp,chosen))

    chosen = chosen[1:,:]
    number_ch = chosen.shape[0]

    plt.plot(chosen[:,0],chosen[:,1],'o', color='black', markersize=0.5)
    plt.title(name+" sig: "+str(sig)+" ave: "+str(ave)+" "+str(number_ch)+' points')
    plt.xlim([0, limits_[0]])
    plt.ylim([0, limits_[1]])
    # plt.show()
    plt.savefig(name+"_"+str(sig)+"_"+str(ave)+'.png')

    print(chosen, chosen.shape)
    print('\n\n')


if __name__ == '__main__':

    sig = 0.1
    ave = 0.5
    # rd = pandas.read_csv('res_'+str(sig)+"_"+str(ave)+'.csv')
    # rd = pandas.read_csv('res_9_21.csv')
    rd = pandas.read_csv('res_9_23.csv')
    num_ex = len(rd['x'])
    dim = int(rd['dim'][1])
    limits_ = list(map(float, re.findall('\d+',(rd['limits_'][1]))))

    extract_plot('chosen', sig, ave)
    extract_plot('data_points', sig, ave)