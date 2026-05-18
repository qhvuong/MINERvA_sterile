import os
import numpy as np

def MergeAsimovs():
    while True:
        print("Type the path of the directory that contains the Asimov delta chi2 files")
        path = input().strip()

        if len(path) == 0:
            print("Type a valid path")
            continue

        if not os.path.isdir(path):
            print("Not a valid path")
            break

        csv_files = [
            f for f in os.listdir(path)
            if f.startswith("sample_dchi2s_") and f.endswith(".csv")
        ]

        if len(csv_files) == 0:
            print("No sample_dchi2s_*.csv files found in {}".format(path))
            break

        # Detect mode from filenames:
        #   sample_dchi2s_noFluxProfile_0.csv
        #   sample_dchi2s_profiledFlux_1.csv
        modes = []
        for f in csv_files:
            tmp = f.replace("sample_dchi2s_", "")
            mode = tmp.rsplit("_", 1)[0]
            modes.append(mode)

        modes = sorted(list(set(modes)))
        print("Found modes:", modes)

        for mode in modes:
            results = []
            files_used = []

            mode_files = [
                f for f in csv_files
                if f.startswith("sample_dchi2s_{}_".format(mode))
            ]

            for f in sorted(mode_files):
                f_path = os.path.join(path, f)
                print("Reading", f_path)

                res = np.loadtxt(f_path, delimiter=",")
                res = np.atleast_1d(res)

                if np.isnan(res).any() or np.isinf(res).any():
                    print("Bad delta chi2 computed in {}".format(f_path))
                    continue

                results.append(res)
                files_used.append(f_path)

            if len(results) == 0:
                print("No valid dchi2 CSV files found for mode:", mode)
                continue

            results = np.concatenate(results)

            print("\nMerged mode:", mode)
            print("  files used:", len(files_used))
            print("  total toys:", results.shape[0])
            print("  min/max/mean:", np.nanmin(results), np.nanmax(results), np.nanmean(results))

            outname = "asimov_deltachi2s_{}.npy".format(mode)
            print("  saving:", outname)
            np.save(outname, results)

        print("Done saving Asimov delta chi2 files")
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

if __name__ == "__main__":
    print("Do you want to merge the chi2 contour files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeChi2s()

    print("Do you want to merge asimov delta chi2 files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeAsimovs()
