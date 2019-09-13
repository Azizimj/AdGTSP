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


def extract_plot(name):
    print(name)
    chosen = np.zeros((1, dim))
    for i, row in enumerate(rd[name]):
        if str(row) != 'nan':
            # print(row)
            tmp = re.findall('\.\d+',row)
            tmp = list(map(float, tmp))
            len_ = len(tmp)
            tmp = np.array(tmp).reshape(int(len_/dim),dim)
            chosen = np.vstack((tmp,chosen))

    chosen = chosen[1:,:]

    plt.plot(chosen[:,0],chosen[:,1],'o', color='black')
    plt.title(name)
    # plt.show()
    plt.savefig(name+'.png')

    print(chosen, chosen.shape)
    print('\n\n')


if __name__ == '__main__':

    rd = pandas.read_csv('res.csv')
    num_ex = len(rd['x'])
    dim = int(rd['dim'][1])

    extract_plot('chosen')
    extract_plot('data_points')