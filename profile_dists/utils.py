import os.path
import psutil
import pandas as pd
import numpy as np
import fastparquet as fp
from numba import jit, njit
from numba.typed import List
import pyarrow.parquet as pq
import pyarrow as pa
import re
import sys
import gc
from profile_dists.constants import MIN_FILE_SIZE, FILE_FORMATS, VALID_INT_TYPES

# Constants
MISSING_ALLELE = '0' # The label for a missing allele
MISSING_ALLELE_DISTANCE = 0 # The distance value for a missing allele

def guess_format(unique_values):
    '''
    Function accepts a list of values and determines if the values are consistent with a list of integers,
    or a set of fixed length hash values
    :param unique_values: list of values
    :return: string (hash, mix, int)
    '''
    length_equal = is_all_same_len(unique_values)
    has_integers = contains_integers(unique_values)
    has_alpha  = contains_alpha(unique_values)

    file_format = ''
    #columns contains hash codes
    if length_equal and has_integers and has_alpha:
        file_format = 'hash'

    # columns contains a mix of integers and other info
    elif has_integers and has_alpha:
        file_format = 'mix'

    #columns contain only integers
    elif has_integers:
        file_format = 'int'

    return file_format


def is_all_same_len(unique_values):
    '''
    Accepts a list of values and calculates the length of each one
    :param unique_values: list of numeric or string values
    :return: True if all values are the same length
    '''
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
    '''
    Accepts a list of values and determines if any represent an integer
    :param unique_values: list of numeric or string values
    :return: True if any value is an integer
    '''
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or re.search('[0-9]+',str(value)):
            status = True
            break
    return status

def contains_alpha(unique_values):
    '''
    Accepts a list of values and determines if any contain [aA -zZ]
    :param unique_values:
    :return:
    '''
    status = False
    for idx, value in enumerate(unique_values):
        if isinstance(value, int) or isinstance(value, float):
            continue
        if re.search('[a-zA-Z]+',str(value)):
            status = True
            break
    return status

def convert_allele_codes(unique_values,method, missing_allele=MISSING_ALLELE, missing_allele_distance=MISSING_ALLELE_DISTANCE):
    '''
    Accepts a list of values and either casts them to an integer for a valid allele or 0 for missing
    :param unique_values: list of numeric or string values
    :param method: [int, hash] to toggle behaviour
    :return: dict of original: integer allele code
    '''
    converted_values = {}
    counter = 1
    for idx,value in enumerate(unique_values):
        if method == 'int':
            if isinstance(value,int) or str(value).isnumeric():
                converted_values[unique_values[idx]] = int(value)
            else:
                converted_values[unique_values[idx]] = missing_allele_distance
        elif method == 'hash':
            if value == missing_allele:
                converted_values[unique_values[idx]] = missing_allele_distance
            else:
                converted_values[unique_values[idx]] = counter
                counter+=1
        else:
            if re.search('[a-zA-Z]+',str(value)) or re.search('[^0-9a-zA-Z]+',str(value)):
                value = missing_allele
            converted_values[unique_values[idx]] = int(value)
    return converted_values


def update_column_map(c1,c2, missing_allele=MISSING_ALLELE, missing_allele_distance=MISSING_ALLELE_DISTANCE):
    '''
    Adds a reference of each column from c2 into c1
    :param c1: dict
    :param c2: dict
    :return: dict
    '''
    if len(c1) == 0:
        allele_id = 1
    else:
        allele_id = max(list(c1.values()))+1

    for k in c2:
        if k == missing_allele:
            c1[k] = missing_allele_distance
        elif not k in c1:
            c1[k] = allele_id
            allele_id+=1

def is_all_columns_int(column_dtypes):
    '''

    :param column_dtypes: List of Pandas column dtypes
    :return: True if all columns are of an integer type
    '''

    for col in column_dtypes:
        if col not in VALID_INT_TYPES:
            return False

    return True

def count_missing_data(df):
    '''
    Scans the df for missing data
    :param df: pandas dataframe
    :return: dict of column: and count of missing values per column
    '''
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
    '''
    Identifies columns with missing data which exceeds the defined threshold
    :param column_counts: dict count of missing data per column
    :param threshold: float or in of the value to filter
    :return: list of columns which fail the test
    '''
    cols_to_remove = []
    for c in column_counts:
        if column_counts[c] > threshold:
            cols_to_remove.append(c)
    return cols_to_remove


def filter_columns(df,columns_to_remove):
    '''
    Removes columns from a data frame
    :param df:
    :param columns_to_remove:
    :return: data frame with columns removed
    '''
    return df.drop(columns_to_remove, axis=1)


def get_header(profile_path,format='text', nrows=10):
    df = pd.DataFrame()
    if format=='text':
        # Note that because of pandas bug, read_csv(index_col=0, dtype=str, ...) will generate
        # an index column that isn't guaranteed to be strings. However, if we're only
        # pulling the headers (which doesn't include the index), then this should be safe.
        df = pd.read_csv(profile_path, header=0, sep="\t", index_col=0, low_memory=False, nrows=nrows, dtype=str)
    elif format=='parquet':
        df = pd.read_parquet(
            profile_path,
            engine='auto',
            columns=None,
            storage_options=None,
            nrows=nrows
        )
    return list(df.columns)

def combine_header(h1,h2):
    combined = h1
    for c in h2:
        if c not in combined:
            combined.append(c)
    return combined

def create_col_map(columns):
    column_mapping = {}
    for col in columns:
        column_mapping[col] = {}
    return column_mapping

def init_combined_header(query_path,ref_path,format='text'):
    h1 = get_header(query_path,format='text', nrows=1)
    h2 = get_header(ref_path,format='text', nrows=1)
    return create_col_map(combine_header(h1,h2))

def process_profile(profile_path,format="text",column_mapping={}, missing_allele=MISSING_ALLELE):
    '''
    Reads in a file in (text, parquet) formats and applies processing to standardize the data and prepare it to be used
    for distance calculations.
    If an allele is missing (NA, ' ', '-', '') then it is replaced with 0
    :param profile_path: path to file
    :param format: format of the file [text, parquet]
    :param column_mapping: Previous allele code mapping to apply to the current file
    :return: (dict, pd)
    '''
    df = pd.DataFrame()
    if format=='text':
        df = pd.read_csv(profile_path, header=0, sep="\t", low_memory=False, dtype=str)

        # There's a bug in some versions of pandas with the read_csv function.
        # If you attempt pd.read_csv(dtype=str, index_col=0),
        # then pandas will not cast the index to the specified type (str).
        # We work around this by not loading an index, and then reshaping the
        # DataFrame to have the index we want, which will be in str format.
        index = df.iloc[:, 0]
        df = df.iloc[:, 1:]
        df = df.set_index(index)

    elif format=='parquet':
        df = pd.read_parquet(
            profile_path,
            engine='auto',
            columns=None,
            storage_options=None,
        )
        # Ensure all values and the index are strings for consistent downstream handling
        df = df.astype(str)
        df.index = df.index.astype(str)
    columns = df.columns.values.tolist()
    if len(column_mapping) > 0:
        missing_fields = list(set(column_mapping.keys()) - set(columns) )
        dtype_change = {}
        for col in missing_fields:
            df[col] = missing_allele
            dtype_change[col] = 'int64'

        header = list(column_mapping.keys())
        df = df.astype(dtype_change)
    else:
        header = columns
    df = df[header]

    column_dtypes = df.dtypes.tolist()
    is_correct_format = is_all_columns_int(column_dtypes)
    #If all columns are already integers then skip the extra processing steps
    if is_correct_format:
        return (column_mapping, df)

    df = df.fillna(missing_allele)
    df = df.replace('?', missing_allele, regex=False)
    df = df.replace(' ', missing_allele, regex=False)
    df = df.replace('-', missing_allele, regex=False)
    df = df.replace('', missing_allele, regex=False)
    df = df.replace('_', missing_allele, regex=False)

    for column in columns:
        unique_col_values = sorted(df[column].unique().tolist())
        method = guess_format(List(unique_col_values))
        converted_allele_codes = convert_allele_codes(unique_col_values, method)
        if not column in column_mapping:
            column_mapping[column] = {}

        update_column_map(column_mapping[column], converted_allele_codes)

        df[column] = df[column].map(column_mapping[column])
    return (column_mapping, df)


def convert_profiles(df):
    '''
    Convert the pd data frame into a numpy array
    :param df: pandas df of integer allele profiles
    :return: (list,list) labels, numpy.arrays
    '''
    labels = [str(x) for x in df.index.tolist()]
    profiles = []
    for index,row in df.iterrows():
        profiles.append(np.array(row.values.tolist()))
    return labels, profiles

def write_profiles(df,out_file,format):
    '''
    Write df to parquet
    :param df:
    :param out_file:
    :param format:
    :return:
    '''
    if format == 'parquet':
        df.to_parquet(out_file,compression='gzip')
    else:
        df.to_csv(out_file,sep="\t",header=True)

@jit(nopython=True)
def count_missing(p):
    '''
    Counts number of 0 occurances in the data
    :param p:
    :return:
    '''
    count = 0
    for idx,value in enumerate(p):
        if value ==0:
            count+=1

    return count


@jit(nopython=True)
def get_distance_raw(p1, p2):
    '''
    Calculates hamming distance with missing data skipped
    :param p1: numpy profile
    :param p2: numpy profile
    :return: int differences
    '''
    count = 0
    for v1,v2 in zip(p1,p2):
        if v1 == 0 or v2 == 0:
            continue
        if v1 != v2:
            count+=1
    return count

@jit(nopython=True)
def get_distance_scaled(p1, p2):
    '''
    Calculates % distance with missing data skipped
    :param p1: numpy profile
    :param p2: numpy profile
    :return: float differences
    '''
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


@jit(nopython=True)
def get_distance_raw_missing(p1, p2):
    '''
    Calculates hamming distance with missing data counted as differences
    :param p1: numpy profile
    :param p2: numpy profile
    :return: int differences
    '''
    count = 0
    for v1, v2 in zip(p1, p2):
        if v1 != v2:
            count += 1
    return count


@jit(nopython=True)
def get_distance_scaled_missing(p1, p2):
    '''
    Calculates % distance with missing data counted as differences
    :param p1: numpy profile
    :param p2: numpy profile
    :return: float differences
    '''
    count_compared_sites = 0
    count_match = 0
    for v1, v2 in zip(p1, p2):
        count_compared_sites += 1
        if v1 == v2:
            count_match += 1

    if count_compared_sites:
        return 100.0 * (float(count_compared_sites) - float(count_match)) / float(count_compared_sites)
    else:
        return 100.0

def calc_batch_size(num_records,num_columns,byte_value_size,max_mem=None):
    '''
    Calculation for the number of simultaneous calcualtions of distance to occur concurrently
    :param num_records: number of records
    :param num_columns: number of columns in a profile
    :param byte_value_size: integer of memory required
    :param max_mem: value to specify the maximum memory to use in the program
    :return: int of the number of records to process concurently
    '''
    if max_mem == None:
        max_mem = psutil.virtual_memory()
        avail = max_mem.available
    else:
        avail = max_mem
    p = (byte_value_size * num_columns) + 56
    profile_mem = p * num_records * 4
    dist_mem = ((byte_value_size * num_records) + 56) ** 2
    estimated_mem_needed = profile_mem + dist_mem
    if estimated_mem_needed < avail:
        return num_records

    avail -= profile_mem
    num_batches = int(avail / ((byte_value_size * num_records) + 56))
    batch_size = int(num_records / num_batches)

    if batch_size <= 0:
        batch_size = 100
    elif batch_size > num_records:
        batch_size = num_records

    return batch_size

@jit(nopython=True)
def validate_file(f):
    '''
    Helper function to determine a file exists and is not empty
    :param f: string path to file
    :return:  True on success
    '''
    if not os.path.isfile(f):
        return False

    if os.path.getsize(f) < MIN_FILE_SIZE:
        return False

    return True

def compare_headers(file1,file2):
    '''
    Compares two files to see if their headers are the same
    :param file1: string path to file 1
    :param file2: string path to file 2
    :return: True on success
    '''
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

def guess_profile_format(f):
    '''
    Helper function to determine what file type a file is
    :param f: string path to file
    :return: string of the format
    '''
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
    '''
    Counts the number of lines in a file
    :param f: string path to file
    :return: int
    '''
    return int(os.popen(f'wc -l {f}').read().split()[0])


def calc_distances_scaled(query_profiles, query_labels, ref_profiles, ref_labels, parquet_file, batch_size=1):
    """
    Calculate pairwise scaled distances between query and reference profiles.

    The distances are expressed as percentages (0–100) and missing data (encoded
    as 0) is ignored when computing the proportion of matching alleles.

    Parameters
    ----------
    query_profiles : list of numpy.ndarray
        List of integer-encoded query profiles (1D arrays of equal length).
    query_labels : list of str
        List of query sample identifiers corresponding to ``query_profiles``.
    ref_profiles : list of numpy.ndarray
        List of integer-encoded reference profiles.
    ref_labels : list of str
        List of reference sample identifiers corresponding to ``ref_profiles``.
    parquet_file : str
        Path to the output parquet file for this batch.
    batch_size : int, optional
        Number of query records to include per I/O batch when writing results.

    Returns
    -------
    None
        Results are written to ``parquet_file`` in parquet format.
    """
    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)

    if num_query_profiles == 0 or num_ref_profiles == 0:
        return

    # Build dense blocks for efficient Numba-based computation.
    q_block = np.vstack(query_profiles).astype(np.int32)
    r_block = np.vstack(ref_profiles).astype(np.int32)

    dist_block = _numba_scaled_block_ignore_missing(q_block, r_block)

    dists = []

    for i in range(num_query_profiles):
        row = [query_labels[i]]
        row.extend(dist_block[i, :].tolist())
        dists.append(row)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression="GZIP")
            else:
                fp.write(parquet_file, df, append=True, compression="GZIP")
            del df
            dists = []
            count = 0
            gc.collect()

    if dists:
        df = pd.DataFrame(dists, columns=columns)
        if not os.path.isfile(parquet_file):
            fp.write(parquet_file, df, compression="GZIP")
        else:
            fp.write(parquet_file, df, append=True, compression="GZIP")

def calc_distances_scaled_missing(query_profiles, query_labels, ref_profiles, ref_labels, parquet_file, batch_size=1):
    """
    Calculate pairwise scaled distances counting missing data as differences.

    The distances are expressed as percentages (0–100) and missing data
    (encoded as 0) is treated the same as any other allele when computing the
    proportion of matching alleles.

    Parameters
    ----------
    query_profiles : list of numpy.ndarray
        List of integer-encoded query profiles (1D arrays of equal length).
    query_labels : list of str
        List of query sample identifiers corresponding to ``query_profiles``.
    ref_profiles : list of numpy.ndarray
        List of integer-encoded reference profiles.
    ref_labels : list of str
        List of reference sample identifiers corresponding to ``ref_profiles``.
    parquet_file : str
        Path to the output parquet file for this batch.
    batch_size : int, optional
        Number of query records to include per I/O batch when writing results.

    Returns
    -------
    None
        Results are written to ``parquet_file`` in parquet format.
    """
    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)

    if num_query_profiles == 0 or num_ref_profiles == 0:
        return

    q_block = np.vstack(query_profiles).astype(np.int32)
    r_block = np.vstack(ref_profiles).astype(np.int32)

    dist_block = _numba_scaled_block_missing_as_diff(q_block, r_block)

    dists = []

    for i in range(num_query_profiles):
        row = [query_labels[i]]
        row.extend(dist_block[i, :].tolist())
        dists.append(row)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression="GZIP")
            else:
                fp.write(parquet_file, df, append=True, compression="GZIP")
            del df
            dists = []
            count = 0
            gc.collect()

    if dists:
        df = pd.DataFrame(dists, columns=columns)
        if not os.path.isfile(parquet_file):
            fp.write(parquet_file, df, compression="GZIP")
        else:
            fp.write(parquet_file, df, append=True, compression="GZIP")

def calc_distances_hamming(query_profiles, query_labels, ref_profiles, ref_labels, parquet_file, batch_size=1):
    """
    Calculate pairwise Hamming distances ignoring missing data.

    Missing data (encoded as 0) is skipped when counting differences between
    alleles.

    Parameters
    ----------
    query_profiles : list of numpy.ndarray
        List of integer-encoded query profiles (1D arrays of equal length).
    query_labels : list of str
        List of query sample identifiers corresponding to ``query_profiles``.
    ref_profiles : list of numpy.ndarray
        List of integer-encoded reference profiles.
    ref_labels : list of str
        List of reference sample identifiers corresponding to ``ref_profiles``.
    parquet_file : str
        Path to the output parquet file for this batch.
    batch_size : int, optional
        Number of query records to include per I/O batch when writing results.

    Returns
    -------
    None
        Results are written to ``parquet_file`` in parquet format.
    """
    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)

    if num_query_profiles == 0 or num_ref_profiles == 0:
        return

    q_block = np.vstack(query_profiles).astype(np.int32)
    r_block = np.vstack(ref_profiles).astype(np.int32)

    dist_block = _numba_hamming_block_ignore_missing(q_block, r_block)

    dists = []

    for i in range(num_query_profiles):
        row = [query_labels[i]]
        row.extend(dist_block[i, :].tolist())
        dists.append(row)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression="GZIP")
            else:
                fp.write(parquet_file, df, append=True, compression="GZIP")
            del df
            dists = []
            count = 0
            gc.collect()

    if dists:
        df = pd.DataFrame(dists, columns=columns)
        if not os.path.isfile(parquet_file):
            fp.write(parquet_file, df, compression="GZIP")
        else:
            fp.write(parquet_file, df, append=True, compression="GZIP")

def calc_distances_hamming_missing(query_profiles, query_labels, ref_profiles, ref_labels, parquet_file, batch_size=1):
    """
    Calculate pairwise Hamming distances counting missing data as differences.

    Missing data (encoded as 0) is treated the same as any other allele when
    counting differences between alleles.

    Parameters
    ----------
    query_profiles : list of numpy.ndarray
        List of integer-encoded query profiles (1D arrays of equal length).
    query_labels : list of str
        List of query sample identifiers corresponding to ``query_profiles``.
    ref_profiles : list of numpy.ndarray
        List of integer-encoded reference profiles.
    ref_labels : list of str
        List of reference sample identifiers corresponding to ``ref_profiles``.
    parquet_file : str
        Path to the output parquet file for this batch.
    batch_size : int, optional
        Number of query records to include per I/O batch when writing results.

    Returns
    -------
    None
        Results are written to ``parquet_file`` in parquet format.
    """
    count = 0
    columns = ["dists"] + [str(x) for x in ref_labels]
    num_query_profiles = len(query_profiles)
    num_ref_profiles = len(ref_profiles)

    if num_query_profiles == 0 or num_ref_profiles == 0:
        return

    q_block = np.vstack(query_profiles).astype(np.int32)
    r_block = np.vstack(ref_profiles).astype(np.int32)

    dist_block = _numba_hamming_block_missing_as_diff(q_block, r_block)

    dists = []

    for i in range(num_query_profiles):
        row = [query_labels[i]]
        row.extend(dist_block[i, :].tolist())
        dists.append(row)
        count += 1

        if count == batch_size:
            df = pd.DataFrame(dists, columns=columns)
            if not os.path.isfile(parquet_file):
                fp.write(parquet_file, df, compression="GZIP")
            else:
                fp.write(parquet_file, df, append=True, compression="GZIP")
            del df
            dists = []
            count = 0
            gc.collect()

    if dists:
        df = pd.DataFrame(dists, columns=columns)
        if not os.path.isfile(parquet_file):
            fp.write(parquet_file, df, compression="GZIP")
        else:
            fp.write(parquet_file, df, append=True, compression="GZIP")

def is_file_ok(f):
    '''
    Helper function to determine if a profile file exists, has a header and >= 1 row of data
    :param f:
    :return: True on success
    '''
    status = True
    if not os.path.isfile(f):
        status = False
    elif get_file_length(f) < 1:
        status = False
    elif os.path.getsize(f) < MIN_FILE_SIZE:
        status = False

    return status

@jit(nopython=True)
def filter_dists(labels,distances,threshold):
    '''
    Filter a list of distances to ones less than or equal to the threshold
    :param labels: list of data labels
    :param distances: list of float/int distances
    :param threshold: float/int
    :return: dict {label: distance} of distances <= threshold
    '''
    results = {}
    for id, value in zip(labels,distances):
        if value <= threshold:
            results[id] = value
    return results


def format_pairwise_dist(df,threshold=-1):
    '''
    Accepts a pd dataframe in matrix format and converts it into query_id,ref_id, distance format df for
    distances <= thresh
    :param df: input df
    :param threshold: float/int
    :return: pd df query_id,ref_id with distances filtered by the threshold
    '''
    dists = {}
    columns = df.columns.values.tolist()[1:]
    for index,row in df.iterrows():
        dists[row[0]] = {row.iloc[0]:0}
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



def write_dist_results(mat, outfile, outtype, outfmt, batch_size=1, threshold=-1):
    """
    Write distance results from the intermediate matrix parquet to the requested
    output format and layout.

    Parameters
    ----------
    mat : str
        Path to the intermediate distance matrix parquet file.
    outfile : str
        Path where the final results file will be written.
    outtype : {"pairwise", "matrix"}
        Logical layout of the output. "matrix" is a wide matrix with one row per
        query and one column per reference. "pairwise" is a long table with
        columns [query_id, ref_id, dist].
    outfmt : {"text", "parquet"}
        File format for the output (TSV text or parquet).
    batch_size : int, optional
        Number of rows to process per parquet batch when streaming.
    threshold : float, optional
        Optional distance threshold; only rows with distances less than or equal
        to this value are retained for pairwise output. Ignored for matrix
        output.
    """
    # If the desired output is a matrix in parquet format simply rename the
    # intermediate matrix file.
    if outtype == 'matrix' and outfmt == 'parquet':
        os.rename(mat, outfile)
        return

    # If the user requested a pairwise parquet file, derive a long-format
    # parquet directly from the intermediate matrix parquet and return.
    if outtype == 'pairwise' and outfmt == 'parquet':
        write_pairwise_parquet_from_matrix_parquet(
            mat,
            outfile,
            threshold=threshold,
            batch_size=batch_size,
        )
        return

    # For all remaining cases we stream over the matrix parquet and emit either
    # a text matrix or a text pairwise file.
    init_file = True
    parquet_file = pq.ParquetFile(mat)

    for batch in parquet_file.iter_batches(batch_size):
        batch_df = batch.to_pandas()
        del batch

        if outtype == 'pairwise':
            batch_df = format_pairwise_dist(batch_df, threshold=threshold)

        if init_file:
            init_file = False
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(outfile, index=False, header=True, sep="\t")
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(outfile, index=False, header=True, sep="\t")
            else:
                # Fallback parquet writer for non-matrix parquet output
                fp.write(outfile, batch_df, compression='GZIP')
        else:
            if outfmt == 'text' and outtype == 'matrix':
                batch_df.to_csv(
                    outfile,
                    mode='a',
                    index=False,
                    header=False,
                    sep="\t",
                )
            elif outfmt == 'text' and outtype == 'pairwise':
                batch_df.to_csv(
                    outfile,
                    mode='a',
                    index=False,
                    header=False,
                    sep="\t",
                )
            else:
                fp.write(outfile, batch_df, append=True, compression='GZIP')


def write_pairwise_parquet_from_matrix_parquet(
    mat, outfile, threshold=-1.0, batch_size=1,
):
    """
    Convert a wide matrix parquet produced by the distance calculation stage
    into a long-format pairwise parquet file.

    The input parquet (``mat``) is expected to have the first column containing
    the query sample identifiers and the remaining columns containing distances
    to each reference sample.

    Parameters
    ----------
    mat : str
        Path to the intermediate distance matrix parquet file.
    outfile : str
        Path where the pairwise parquet file will be written.
    threshold : float, optional
        Optional distance threshold; only distances less than or equal to this
        value are written when a non-negative value is provided. When -1
        (default), all distances are reported.
    batch_size : int, optional
        Number of rows to process per parquet batch when streaming.
    """
    parquet_file = pq.ParquetFile(mat)

    schema = pa.schema(
        [
            ("query_id", pa.string()),
            ("ref_id", pa.string()),
            ("dist", pa.float32()),
        ]
    )

    writer = None
    try:
        for batch in parquet_file.iter_batches(batch_size):
            columns = batch.schema.names
            if not columns:
                continue

            # First column contains the query/sample identifiers
            query_ids = batch.column(0).to_numpy(zero_copy_only=False)

            # Remaining columns correspond to reference IDs
            ref_ids = np.array(columns[1:], dtype=object)
            if ref_ids.size == 0:
                continue

            # Build a dense 2D NumPy array of distances with shape (rows, n_refs)
            dist_cols = []
            for col_idx in range(1, len(columns)):
                col = batch.column(col_idx).to_numpy(zero_copy_only=False)
                dist_cols.append(col)
            dist_block = np.vstack(dist_cols).T

            q_list = []
            r_list = []
            d_list = []

            for i, qid in enumerate(query_ids):
                row = dist_block[i, :]
                if threshold >= 0:
                    mask = row <= threshold
                else:
                    mask = np.ones_like(row, dtype=bool)

                if not mask.any():
                    continue

                selected_refs = ref_ids[mask]
                selected_dists = row[mask]

                q_list.extend([str(qid)] * len(selected_refs))
                r_list.extend(selected_refs.astype(str).tolist())
                d_list.extend(selected_dists.astype("float32").tolist())

            if not q_list:
                continue

            table = pa.table(
                {
                    "query_id": pa.array(q_list, type=pa.string()),
                    "ref_id": pa.array(r_list, type=pa.string()),
                    "dist": pa.array(d_list, type=pa.float32()),
                },
                schema=schema,
            )

            if writer is None:
                writer = pq.ParquetWriter(outfile, schema=schema)
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
def get_missing_loci_counts(profiles,labels,count_loci):
    '''
    Counts the % missing data by sample
    :param profiles: list numpy arrays
    :param labels: list labels
    :param count_loci: total number of loci in profile
    :return: dict {label: % missing}
    '''
    n = len(labels)
    counts = {}
    for i in range(0, n):
        counts[labels[i]] = count_missing(profiles[i])/count_loci
    return counts

def flag_samples(missing_counts,threshold):
    '''
    Produces a list of labels which exceed the threshold for missing data
    :param missing_counts: dict of {label: % missing}
    :param threshold: int/float
    :return: list of labels
    '''
    r = []
    for sample_id in missing_counts:
        if missing_counts[sample_id] > threshold:
            r.append(sample_id)
    return sorted(r)

def filter_samples(labels,profiles,labels_to_remove):
    '''
    Removes samples from profile which are flagged to be removed
    :param labels: list of labels
    :param profiles: list of numpy profiles
    :param labels_to_remove:  List of labels to remove
    :return: (list, list) labels, list numpy profiles
    '''
    l = []
    p = []
    for idx,label in enumerate(labels):
        if label in labels_to_remove:
            continue
        l.append(label)
        p.append(profiles[idx])
    return l, p



@njit
def _numba_scaled_block_ignore_missing(q_block, r_block):
    """Compute scaled percent distance (0–100) ignoring missing (value 0).

    Parameters
    ----------
    q_block : numpy.ndarray
        2D array of integer-encoded query profiles with shape (n_query, n_loci).
    r_block : numpy.ndarray
        2D array of integer-encoded reference profiles with shape (n_ref, n_loci).

    Returns
    -------
    numpy.ndarray
        2D array of float32 distances with shape (n_query, n_ref).
    """
    n_query, n_loci = q_block.shape
    n_ref = r_block.shape[0]
    out = np.empty((n_query, n_ref), dtype=np.float32)
    for i in range(n_query):
        for j in range(n_ref):
            count_compared = 0
            count_match = 0
            for k in range(n_loci):
                v1 = q_block[i, k]
                v2 = r_block[j, k]
                if v1 == 0 or v2 == 0:
                    continue
                count_compared += 1
                if v1 == v2:
                    count_match += 1
            if count_compared > 0:
                out[i, j] = 100.0 - (100.0 * (count_match / count_compared))
            else:
                out[i, j] = 100.0
    return out


@njit
def _numba_hamming_block_ignore_missing(q_block, r_block):
    """Compute Hamming distances ignoring missing (value 0).

    Parameters
    ----------
    q_block : numpy.ndarray
        2D array of integer-encoded query profiles with shape (n_query, n_loci).
    r_block : numpy.ndarray
        2D array of integer-encoded reference profiles with shape (n_ref, n_loci).

    Returns
    -------
    numpy.ndarray
        2D array of float32 distances with shape (n_query, n_ref).
    """
    n_query, n_loci = q_block.shape
    n_ref = r_block.shape[0]
    out = np.empty((n_query, n_ref), dtype=np.float32)
    for i in range(n_query):
        for j in range(n_ref):
            diff_count = 0
            for k in range(n_loci):
                v1 = q_block[i, k]
                v2 = r_block[j, k]
                if v1 == 0 or v2 == 0:
                    continue
                if v1 != v2:
                    diff_count += 1
            out[i, j] = diff_count
    return out


@njit
def _numba_scaled_block_missing_as_diff(q_block, r_block):
    """Compute scaled percent distance (0–100) counting missing as differences.

    Parameters
    ----------
    q_block : numpy.ndarray
        2D array of integer-encoded query profiles with shape (n_query, n_loci).
    r_block : numpy.ndarray
        2D array of integer-encoded reference profiles with shape (n_ref, n_loci).

    Returns
    -------
    numpy.ndarray
        2D array of float32 distances with shape (n_query, n_ref).
    """
    n_query, n_loci = q_block.shape
    n_ref = r_block.shape[0]
    out = np.empty((n_query, n_ref), dtype=np.float32)
    for i in range(n_query):
        for j in range(n_ref):
            count_compared = 0
            count_match = 0
            for k in range(n_loci):
                v1 = q_block[i, k]
                v2 = r_block[j, k]
                count_compared += 1
                if v1 == v2:
                    count_match += 1
            if count_compared > 0:
                out[i, j] = 100.0 - (100.0 * (count_match / count_compared))
            else:
                out[i, j] = 100.0
    return out


@njit
def _numba_hamming_block_missing_as_diff(q_block, r_block):
    """Compute Hamming distances counting missing as differences.

    Parameters
    ----------
    q_block : numpy.ndarray
        2D array of integer-encoded query profiles with shape (n_query, n_loci).
    r_block : numpy.ndarray
        2D array of integer-encoded reference profiles with shape (n_ref, n_loci).

    Returns
    -------
    numpy.ndarray
        2D array of float32 distances with shape (n_query, n_ref).
    """
    n_query, n_loci = q_block.shape
    n_ref = r_block.shape[0]
    out = np.empty((n_query, n_ref), dtype=np.float32)
    for i in range(n_query):
        for j in range(n_ref):
            diff_count = 0
            for k in range(n_loci):
                v1 = q_block[i, k]
                v2 = r_block[j, k]
                if v1 != v2:
                    diff_count += 1
            out[i, j] = diff_count
    return out

