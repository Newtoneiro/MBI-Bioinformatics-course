import argparse
import pandas as pd

def main(path_finalcall, path_dgv):
    dgv_data = pd.read_csv(path_dgv, header=0, sep='\t',
                         dtype={'chr': str, 'start': int, 'end': int, 'variantsubtype': str, 'samples': str})
    dgv_data = dgv_data[['chr', 'start', 'end', 'variantsubtype', 'samples']]
    dgv_data = dgv_data[dgv_data['chr'] == '20']

    cnv_data = pd.read_csv(path_finalcall, sep=',')
    cnv_data = cnv_data[['sample_name', 'chr', 'st_bp', 'ed_bp', 'cnv']]

    deletion_overlapping_count, duplication_overlapping_count, overlapping_more_than_80_count = 0, 0, 0

    for _, cnv_row in cnv_data.iterrows():
        print(cnv_row['sample_name']) # show sample name
        x1, x2 = cnv_row['st_bp'], cnv_row['ed_bp']
        for _, dgv_row in dgv_data.iterrows():
            y1, y2 = dgv_row['start'], dgv_row['end']
            if max(x1, y1) <= min(x2, y2):
                if dgv_row['variantsubtype'] == 'deletion':
                    deletion_overlapping_count += 1
                if dgv_row['variantsubtype'] == 'duplication':
                    duplication_overlapping_count += 1
                if min(x2, y2) - max(x1, y1) >= 0.8 * (y2 - y1):
                    overlapping_more_than_80_count += 1

    print(f'Overlapping duplications: {duplication_overlapping_count}')
    print(f'Overlapping deletions: {deletion_overlapping_count}')
    print(f'Any changes overlapping by more than 80%: {overlapping_more_than_80_count}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('finalcall', type=str, help='Path to finalcall file')
    parser.add_argument('dgv', type=str, help='Path to dgv file')

    args = parser.parse_args()

    main(args.finalcall, args.dgv)
