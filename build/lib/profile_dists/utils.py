import os.path
import shutil
import sys
import time
import psutil
import pandas as pd
import numpy as np
import fastparquet as fp
import tables
from numba import jit
from numba.typed import List
import pyarrow.parquet as pq
import re
from profile_dists.constants import MIN_FILE_SIZE, FILE_FORMATS, VALID_INT_TYPES


def guess_format(unique_values):
    length_equal = is_all_same_len(unique_values)
    has_integers = contains_integers(unique_values)
    has_alpha  = contains_alpha(unique_values)

    format = ''
    #columns contains hash codes
    if length_equal and has_integers and has_alpha:
        format = 'hash'

    # columns contains a mix of integers and other info
    elif has_integers and has_alpha:
        format = 'mix'

    #columns contain only integers
    elif has_integers:
        format = 'int'

    return format


def is_all_same_len(unique_values):
    l = set()
    for idx,value in enumerate(unique_values):
        if value != '0':
            l.add(len(str(value)))
    if len(l) == 1:
        status = True
    else:
        status = False
    return status


def contains_integers(unique_values):
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or re.search('[0-9]+',value):
            status = True
            break
    return status

def contains_alpha(unique_values):
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or isinstance(value, float):
            continue
        if re.search('[a-zA-Z]+',value):
            status = True
            break
    return status

def convert_allele_codes(unique_values,method):
    converted_values = {}
    counter = 1
    for idx,value in enumerate(unique_values):
        if method == 'int':
            converted_values[unique_values[idx]] = int(value)
        elif method == 'hash':
            if value == '0':
                converted_values[unique_values[idx]] = 0
            else:
                converted_values[unique_values[idx]] = counter
                counter+=1
        else:
            if re.search('[a-zA-Z]+',value) or re.search('\.|~|-',value):
                value = '0'
            converted_values[unique_values[idx]] = int(value)
    return converted_values


def update_column_map(c1,c2):
    for k in c2:
        if not k in c1:
            c1[k] = c2[k]

def is_all_columns_int(column_dtypes):
    count_non_int = 0
    for col in column_dtypes:
        if col in VALID_INT_TYPES:
            continue
        count_non_int+=1
    if count_non_int > 0:
        return False
    return True

def count_missing_data(df):
    counts = {}
    columns = df.columns.values.tolist()
    for c in columns:
        counts[c] = 0
        v = df[c].value_counts()
        if 0 in v:
            counts[c] = v[0]
    return counts

@jit(nopython=True)
def identify_cols_to_remove(column_counts,threshold):
    cols_to_remove = []
    for c in column_counts:
        if column_counts[c] > threshold:
            cols_to_remove.append(c)
    return cols_to_remove


def filter_columns(df,columns_to_remove):
    return df.drop(columns_to_remove, axis=1)


def process_profile(profile_path,format="text",column_mapping={}):

    if format=='text':
        df = pd.read_csv(profile_path,header=0,sep="\t",index_col=0,low_memory=False)
    elif format=='parquet':
        df = pd.read_parquet(
            profile_path,
            engine='auto',
            columns=None,
            storage_options=None,
        )

    columns = df.columns.values.tolist()
    column_dtypes = df.dtypes.tolist()
    is_correct_format = is_all_columns_int(column_dtypes)

    #If all columns are already integers then skip the extra processing steps
    if is_correct_format:
        return (column_mapping, df)

    df = df.replace('?', '0', regex=False)
    df = df.replace(' ', '0', regex=False)
    df = df.replace('-', '0', regex=False)
    df = df.replace('', '0', regex=False)

    for column in columns:
        unique_col_values = sorted(df[column].unique().tolist())
        method = guess_format(List(unique_col_values))
        if not column in column_mapping:
            column_mapping[column] = convert_allele_codes(unique_col_values, method)
        else:
            update_column_map(column_mapping[column], convert_allele_codes(unique_col_values, method))

        df[column] = df[column].map(column_mapping[column])
    return (column_mapping, df)


def convert_profiles(df):
    labels = df.index.tolist()
    profiles = []
    for index,row in df.iterrows():
        profiles.append(np.array(row.values.tolist()))
    return labels, profiles

def write_profiles(df,out_file,format):
    if format == 'parquet':
        df.to_parquet(out_file,compression='gzip')
    else:
        df.to_csv(out_file,sep="\t",header=True)

@jit(nopython=True)
def count_missing(p):
    count = 0
    for idx,value in enumerate(p):
        if value ==0:
            count+=1

    return count


@jit(nopython=True)
def get_distance_raw(p1, p2):
    count = 0
    for v1,v2 in zip(p1,p2):
        if v1 == 0 or v2 == 0:
            continue
        if v1 != v2:
            count+=1
    return count

@jit(nopython=True)
def get_distance_scaled(p1, p2):
    count_compared_sites = 0
    count_match = 0
    for v1,v2 in zip(p1,p2):
        if v1 == 0 or v2 == 0:
            continue
        count_compared_sites+=1
        if v1 == v2:
            count_match+=1
    if count_compared_sites:
        return 100.0 * (float(count_compared_sites) - float(count_match)) / float(count_compared_sites)
    else:
        return 100.0


def calc_batch_size(num_records,num_columns,byte_value_size):
    mem = psutil.virtual_memory()
    avail = mem.available
    p = (byte_value_size * num_columns) + 56
    estimated_mem_needed = p * num_records
    if estimated_mem_needed < avail:
        return num_records
    return int(avail / p)

@jit(nopython=True)
def validate_file(f):
    if not os.path.isfile(f):
        return False

    if os.path.getsize(f) < MIN_FILE_SIZE:
        return False

    return True

def compare_headers(file1,file2):
    h1 = []
    h2 = []
    with open(file1,'r') as f1:
        h1 = next(f1).rstrip().split("\t")
        with open(file2, 'r') as f2:
            h2 = next(f2).rstrip().split("\t")
    if len(h1) > 0 and len(h2) > 0 and len(h1) == len(h2):
        ovl = set(h1) & set(h2)
        if len(ovl) == len(h1):
            return True
    return False

@jit(nopython=True)
def guess_profile_format(f):
    ext = FILE_FORMATS
    ftype = ''

    for format in ext:
        for e in ext[format]:
            if f.endswith(e):
                ftype = format
                break
        if ftype != '':
            break

    return ftype


def get_file_length(f):
    return int(os.popen(f'wc -l {f}').read().split()[0])


def calc_distances_scaled(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):

    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            d.append(get_distance_scaled(query_profiles[i], ref_profiles[k]))
        dists.append(d)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression='GZIP')
            else:
                fp.write(parquet_file, df, append=True, compression='GZIP')
            dists = []
            count = 0

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')

def calc_distances_hamming(query_profiles,query_labels,ref_profiles,ref_labels,parquet_file,batch_size=1):
    count = 0
    columns = ["dists"] + ref_labels
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)
    dists = []

    #Clear an existing file as this can cause unexpected behaviour
    if os.path.isfile(parquet_file):
        os.remove(parquet_file)

    for i in range(0, num_query_profiles):
        d = [ query_labels[i] ]
        for k in range(0, num_ref_profiles):
            d.append(get_distance_raw(query_profiles[i], ref_profiles[k]))
        dists.append(d)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression='GZIP')
            else:
                fp.write(parquet_file, df, append=True, compression='GZIP')
            dists = []
            count = 0

    df = pd.DataFrame(dists, columns=columns)
    if not os.path.isfile(parquet_file):
        fp.write(parquet_file, df, compression='GZIP')
    else:
        fp.write(parquet_file, df, append=True, compression='GZIP')


def is_file_ok(f):
    status = True
    if not os.path.isfile(f):
        status = False
    elif get_file_length(f) < 2:
        status = False
    elif os.path.getsize(f) < MIN_FILE_SIZE:
        status = False

    return status

@jit(nopython=True)
def filter_dists(labels,distances,threshold):
    results = {}
    for id, value in zip(labels,distances):
        if value <= threshold:
            results[id] = value
    return results


def format_pairwise_dist(df,threshold=-1):
    dists = {}
    columns = df.columns.values.tolist()[1:]
    for index,row in df.iterrows():
        dists[row[0]] = {row[0]:0}
        if threshold != -1:
            dists[row[0]] = filter_dists(List(columns), List(row[1:]), threshold)
        else:
            dists[row[0]] = dict(zip(columns, row[1:]))
        dists[row[0]] = {k: v for k, v in sorted(dists[row[0]].items(), key=lambda item: item[1])}

    results = {
        'query_id':[],
        'ref_id':[],
        'dist':[]
    }

    for qid in dists:
        results['query_id'] += [qid] * len(dists[qid])
        results['ref_id'] += list(dists[qid].keys())
        results['dist'] += list(dists[qid].values())


    return pd.DataFrame(results)


def write_dist_results(mat,outfile,outtype,outfmt,batch_size=1,threshold=-1):

    #If the desired output is a matrix in parquet format simply rename the mat file
    if outtype == 'matrix' and outfmt == 'parquet':
        os.rename(mat,outfile)
        return
    init_file = True
    parquet_file = pq.ParquetFile(mat)
    for batch in parquet_file.iter_batches(batch_size):
        batch_df = batch.to_pandas()

        if outtype == 'pairwise':
            batch_df = format_pairwise_dist(batch_df, threshold=threshold)
        if init_file:
            init_file = False
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(outfile,index = False, header = True, sep="\t")
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(outfile, index=False, header=True, sep="\t")
            else:
                if not os.path.isfile(outfile):
                    fp.write(outfile, batch_df, compression='GZIP')
        else:
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(outfile, mode ='a', index = False, header = False, sep="\t")
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(outfile, mode ='a', index = False, header = False, sep="\t")
            else:
                fp.write(parquet_file, batch_df, append=True, compression='GZIP')

def get_missing_loci_counts(profiles,labels):
    n = len(labels)
    counts = {}
    for i in range(0, n):
        counts[labels[i]] = count_missing(profiles[i])/n
    return counts

def flag_samples(missing_counts,threshold):
    r = []
    for sample_id in missing_counts:
        if missing_counts[sample_id] > threshold:
            r.append(sample_id)
    return sorted(r)

def filter_samples(labels,profiles,labels_to_remove):
    l = []
    p = []
    for idx,label in enumerate(labels):
        if label in labels_to_remove:
            continue
        l.append(label)
        p.append(profiles[idx])
    return l, p


#write_dist_results('/Users/jrobertson/Desktop/ListeraClusterComp/BioNumerics.alleles.profiles5.parquet','/Users/jrobertson/Desktop/ListeraClusterComp/BioNumerics.alleles.profiles5.csv','pairwise','text',batch_size=500,threshold=10)