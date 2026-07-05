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
        path = input().strip()

        if len(path) == 0:
            print("Enter a valid path")
            break

        if path[-1] != '/':
            path = path + '/'

        if not os.path.isdir(path):
            print("directory does not exist")
            break

        data_files = [
            f for f in os.listdir(path)
            if f.startswith("chi2_surface_data_") and f.endswith(".npy")
        ]

        if len(data_files) == 0:
            print("No chi2 surface files found in {}".format(path))
            break

        # Detect mode from filenames:
        # chi2_surface_data_noFluxProfile_m_0.1.npy
        modes = []
        for f in data_files:
            tmp = f.replace("chi2_surface_data_", "")
            mode = tmp.split("_m_")[0]
            modes.append(mode)

        modes = sorted(list(set(modes)))
        print("Found modes:", modes)

        for mode in modes:
            mode_data_files = [
                f for f in data_files
                if f.startswith("chi2_surface_data_{}_m_".format(mode))
            ]

            m_names = []
            for f in mode_data_files:
                m_str = f.replace("chi2_surface_data_{}_m_".format(mode), "")
                m_str = m_str.replace(".npy", "")
                m_names.append(float(m_str))

            m_names = sorted(list(set(m_names)))

            datas = []
            asimovs = []
            datas_penalty = []
            asimovs_penalty = []
            good_m_names = []

            for m in m_names:
                data_file = path + "chi2_surface_data_{}_m_{}.npy".format(mode, m)
                data_penalty_file = path + "chi2_penalty_data_{}_m_{}.npy".format(mode, m)

                asimov_file = path + "chi2_surface_pseudodata_{}_m_{}.npy".format(mode, m)
                asimov_penalty_file = path + "chi2_penalty_pseudodata_{}_m_{}.npy".format(mode, m)

                files = [
                    data_file,
                    data_penalty_file,
                    asimov_file,
                    asimov_penalty_file,
                ]

                missing = [f for f in files if not os.path.isfile(f)]
                if len(missing) > 0:
                    print("Missing files for m = {}".format(m))
                    for f in missing:
                        print("  missing:", f)
                    continue

                data = np.load(data_file)
                data_penalty = np.load(data_penalty_file)
                asimov = np.load(asimov_file)
                asimov_penalty = np.load(asimov_penalty_file)

                bad = (
                    np.isnan(data).any() or np.isinf(data).any() or
                    np.isnan(asimov).any() or np.isinf(asimov).any() or
                    np.isnan(data_penalty).any() or np.isinf(data_penalty).any() or
                    np.isnan(asimov_penalty).any() or np.isinf(asimov_penalty).any()
                )

                if bad:
                    print("Bad values found for m = {}".format(m))
                    continue

                datas.append(data)
                datas_penalty.append(data_penalty)

                asimovs.append(asimov)
                asimovs_penalty.append(asimov_penalty)

                good_m_names.append(m)

            datas = np.array(datas)
            datas_penalty = np.array(datas_penalty)

            asimovs = np.array(asimovs)
            asimovs_penalty = np.array(asimovs_penalty)

            good_m_names = np.array(good_m_names)

            print("Merged mode:", mode)
            print("  data shape          =", datas.shape)
            print("  data penalty shape  =", datas_penalty.shape)
            print("  asimov shape        =", asimovs.shape)
            print("  asimov penalty shape=", asimovs_penalty.shape)
            print("  mass slices         =", good_m_names.shape[0])

            np.save("data_chi2s_{}".format(mode), datas)
            np.save("data_penalties_{}".format(mode), datas_penalty)

            np.save("asimov_chi2s_{}".format(mode), asimovs)
            np.save("asimov_penalties_{}".format(mode), asimovs_penalty)

            np.save("delta_m_values_{}".format(mode), good_m_names)

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
