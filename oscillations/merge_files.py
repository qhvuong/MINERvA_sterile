import os
import numpy as np

def MergeAsimovs():
    while True:
        print("Type the path of the directory that contains the Asimov delta chi2 files")
        path = input()
        if len(path) == 0:
            print("Type a valid path")
            continue

        if os.path.isdir(path):
            results = []
            for f in os.listdir(path):
                f_path = path+'/'+f
                if os.path.isdir(f_path):
                    for f_ in os.listdir(f_path):
                        file = f_path+'/'+f_
                        res = np.loadtxt(file,delimiter=',')
                        if np.isnan(res).any() or np.isinf(res).any():
                            print("Bad delta chi2 computed in {}".format(file))
                            continue
                        results.append(res)
                else:
                    res = np.loadtxt(f_path,delimiter=',')
                    if np.isnan(res).any() or np.isinf(res).any():
                        print("Bad delta chi2 computed in {}".format(f_path))
                        continue
                    results.append(res)
            results = np.array(results).flatten()
            print('saving asimov_deltachi2s.npy with {} entries'.format(results.shape[0]))
            np.save("asimov_deltachi2s",results)
            break
        else:
            print("Not a valid path")
            break

def MergeChi2s():
    while True:
        print("Type the path of the directory that contains the chi2 contours you want to merge")
        path = input()
        if len(path) == 0:
            print("Enter a valid path")
            break

        if path[-1] != '/':
            path = path+'/'
            
        if os.path.isdir(path):
            data_files = [f for f in os.listdir(path) if "chi2_surface" in f]
            if len(data_files) == 0:
                print("No files found in {}".format(path))
                break
        else:
            print("directory does not exist")
            break

        # ----- sort file names to order PMNS parametesr ----- #
        start = '_m_'
        end = '_Ue4'
        m_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        m_names = list(set(m_names))
        m_names.sort()

        start = 'Ue4_'
        end = '.dat.npy'
        e_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        e_names = list(set(e_names))
        e_names.sort()

        datas = []
        asimovs = []
        datas_penalty = []
        asimovs_penalty = []

        for m in m_names:
            m_data   = []
            m_data_penalty = []

            m_asimov = []
            m_asimov_penalty = []

            for e in e_names:
                data = np.load(path+"chi2_surface_data_m_{}_Ue4_{}.dat.npy".format(m,e))
                data_penalty = np.load(path+"chi2_penalty_data_m_{}_Ue4_{}.dat.npy".format(m,e))

                asimov = np.load(path+"chi2_surface_pseudodata_m_{}_Ue4_{}.dat.npy".format(m,e))
                asimov_penalty = np.load(path+"chi2_penalty_pseudodata_m_{}_Ue4_{}.dat.npy".format(m,e))

                m_data.append(data)
                m_data_penalty.append(data_penalty)

                m_asimov.append(asimov)
                m_asimov_penalty.append(asimov_penalty)

            datas.append(m_data)
            datas_penalty.append(m_data_penalty)

            asimovs.append(m_asimov)
            asimovs_penalty.append(m_asimov_penalty)

        datas = np.array(datas)
        datas_penalty = np.array(datas_penalty)

        asimovs = np.array(asimovs)
        asimovs_penalty = np.array(asimovs_penalty)

        np.save("data_chi2s",datas)
        np.save("data_penalties",datas_penalty)
        
        np.save("asimov_chi2s",asimovs)
        np.save("asimov_penalties",asimovs_penalty)

        print("Done saving files")

        break

if __name__ in "__main__":
    print("Do you want to merge the chi2 contour files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeChi2s()

    print("Do you want to merge asimov delta chi2 files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeAsimovs()
