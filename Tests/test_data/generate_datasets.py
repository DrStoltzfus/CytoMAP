import numpy as np
import pandas as pd
import os


def gen_label(num, begin='A'):
    return '{} {}'.format(begin, num)


def generate_group(num_common=15, num_channels=20, num_csvs=5, is_3d=True, begin='A', save_dir='Data3D/dataset_a'):
    num_generate_channels = max(0, num_channels - num_common)
    common_channels = [gen_label(idx, begin) for idx in range(num_common - 3)]
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    uncommon_idx = num_common
    for idx in range(num_csvs):
        data = np.random.rand(1000, num_channels)
        data[:, :5] = data[:, :5] * 1e3
        data[:, 5:10] = data[:, 5:10] * 1e1
        data[:, 10:15] = data[:, 10:15] * 1e-3
        data[:, 15:20] = data[:, 15:20] * 1e5
        if not is_3d:
            data[:, 2] = 0

        data = pd.DataFrame(data=data, columns=['X', 'Y', 'Z'] + common_channels +
                                               [gen_label(jdx, begin) for jdx in range(
                                                   uncommon_idx, uncommon_idx + num_generate_channels)
                                                ]
                            )

        uncommon_idx += num_generate_channels
        data.to_csv('{}/{}.csv'.format(save_dir, idx))


generate_group()
generate_group(is_3d=False, save_dir='Data2D/dataset_a', begin='B', num_channels=22)
