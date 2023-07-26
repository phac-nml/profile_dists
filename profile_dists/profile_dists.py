import copy
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import json
import os
from datetime import datetime
from profile_dists.version import __version__
from profile_dists.utils import process_profile, is_file_ok, filter_columns, \
    count_missing_data, write_profiles, convert_profiles, calc_distances_scaled, calc_distances_hamming, \
    calc_distances_scaled_missing, calc_distances_hamming_missing,\
    write_dist_results, calc_batch_size, get_missing_loci_counts, flag_samples, filter_samples

from profile_dists.constants import RUN_DATA

def parse_args():
    """ Argument Parsing method.

        A function to parse the command line arguments passed at initialization of Clade-o-matic,
        format these arguments,  and return help prompts to the user shell when specified.

        Returns
        -------
        ArgumentParser object
            The arguments and their user specifications, the usage help prompts and the correct formatting
            for the incoming argument (str, int, etc.)
        """
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        """
                Class to instantiate the formatter classes required for the argument parser.
                Required for the correct formatting of the default parser values

                Parameters
                ----------
                ArgumentDefaultsHelpFormatter object
                    Instatiates the default values for the ArgumentParser for display on the command line.
                RawDescriptionHelpFormatter object
                    Ensures the correct display of the default values for the ArgumentParser
                """
        pass

    parser = ArgumentParser(
        description="Profile Dists: Calculate genetic distances based on allele profiles v. {}".format(__version__),
        formatter_class=CustomFormatter)
    parser.add_argument('--query','-q', type=str, required=True, help='Query allelic profiles')
    parser.add_argument('--ref','-r', type=str, required=True, help='Reference allelic profiles')
    parser.add_argument('--outdir', '-o', type=str, required=True, help='Result output files')
    parser.add_argument('--outfmt', '-u', type=str, required=False, help='Out format [matrix, pairwise]',default='matrix')
    parser.add_argument('--file_type', '-e', type=str, required=False, help='Out format [text, parquet]',default='text')
    parser.add_argument('--distm', '-d', type=str, required=False, help='Distance method raw hamming or scaled difference [hamming, scaled]',default='scaled')
    parser.add_argument('--missing_thresh', '-t', type=float, required=False,
                        help='Maximum percentage of missing data allowed per locus (0 - 1)',default=1.0)
    parser.add_argument('--sample_qual_thresh', '-c', type=float, required=False,
                        help='Maximum percentage of missing data allowed per sample (0 - 1)',default=1.0)
    parser.add_argument('--match_threshold', '-a', type=str, required=False,
                        help='Either a integer or float depending on what distance method is used (only used with pairwise format',default=-1)
    parser.add_argument('--mapping_file', '-m', type=str, required=False,
                        help='json formatted allele mapping')
    parser.add_argument('--batch_size',  type=int, required=False,
                        help='Manual selection of how many records should be included in a batch, default=auto')
    parser.add_argument('--max_mem',  type=int, required=False,
                        help='Maximum amount of memory to use',default=None)
    parser.add_argument('--force','-f', required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument('-s', '--skip', required=False, help='Skip QA/QC steps',
                        action='store_true')
    parser.add_argument('-n', '--count_missing', required=False, help='Count missing as differences',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()


def main():
    cmd_args = parse_args()
    query_profile = cmd_args.query
    ref_profile = cmd_args.ref
    outdir = cmd_args.outdir
    outfmt = cmd_args.outfmt
    file_type = cmd_args.file_type
    dist_method = cmd_args.distm
    missing_threshold = cmd_args.missing_thresh
    allele_mapping_file = cmd_args.mapping_file
    force = cmd_args.force
    match_threshold = cmd_args.match_threshold
    sample_qual_thresh = cmd_args.sample_qual_thresh
    skip = cmd_args.skip
    count_missing_sites = cmd_args.count_missing
    max_mem = cmd_args.max_mem
    batch_size = cmd_args.batch_size

    run_data = RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    input_files = [query_profile,ref_profile,allele_mapping_file]
    for f in input_files:
        if f is None:
            continue
        if not is_file_ok(f):
            print(f'file {f} either does not exist or is too small to be valid')
            sys.exit()

    allele_map = {}
    if allele_mapping_file is not None:
        with open(allele_mapping_file) as mapping_fh:
            allele_map = json.loads(mapping_fh.read())

    if not force and os.path.isdir(outdir):
        print(f'folder {outdir} already exists, please choose new directory or use --force')
        sys.exit()

    if os.path.isdir(outdir):
        print(f'folder {outdir} already exists, and force specified, cleaning up directory')
        files = [os.path.join(outdir,allele_map.json),
                 os.path.join(outdir, f'query_profile.{file_type}'),
                 os.path.join(outdir, f'ref_profile.{file_type}'),
                 os.path.join(outdir,f'dists.parquet'),
                 os.path.join(outdir, f'results.{file_type}'),
                 os.path.join(outdir, 'run.json'),]
        for file in files:
            if os.path.isfile(file):
                os.remove(file)


    if outfmt != 'matrix' and outfmt != 'pairwise':
        print(f'Supplied format does not match [matrix,pairwise]: {outfmt} ')
        sys.exit()

    if not file_type in ['text', 'parquet']:
        print(f'Supplied filetype does not match [text, parquet]: {outfmt} ')
        sys.exit()

    if not dist_method  in ['hamming','scaled']:
        print(f'Supplied filetype does not match [hamming, scaled]: {dist_method} ')
        sys.exit()

    if missing_threshold < 0 or missing_threshold > 1:
        print(f'Supplied threshold is not between 0 - 1: {missing_threshold} ')
        sys.exit()

    # initialize analysis directory
    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    print(f'Reading query profile: {query_profile}')
    sys.stdout.flush()
    (allele_map, qdf) = process_profile(query_profile,column_mapping=allele_map)


    #Skip calculation for symmetic matrix calculation on same set of data
    if query_profile == ref_profile:
        print(f'Reusing data from: {query_profile}')
        rdf = copy.deepcopy(qdf)
    else:
        print(f'Reading reference profile: {query_profile}')
        (allele_map, rdf) = process_profile(ref_profile, column_mapping=allele_map)

    print(f'Writting allele map')
    sys.stdout.flush()
    with open(os.path.join(outdir,"allele_map.json"),'w' ) as fh:
        fh.write(json.dumps(allele_map, indent=4))

    qcols = set(qdf.columns.values.tolist())
    rcols = set(rdf.columns.values.tolist())
    common_cols = sorted(list(qcols & rcols))

    if len(common_cols) == 0:
        print(f'Error there are no columns in common between: {query_profile}\t{ref_profile}')
        sys.exit()

    #remove cols not present in both
    qcols_to_remove = qcols - set(common_cols)
    run_data['loci_removed'] = list(qcols_to_remove)

    if len(qcols_to_remove) > 0:
        qdf = filter_columns(qdf, qcols_to_remove)

    rcols_to_remove = rcols - set(common_cols)
    if len(rcols_to_remove) > 0:
        rdf = filter_columns(rdf, qcols_to_remove)

    cols_to_remove = []
    if not skip:
        print(f'Performing QA/QC on input allele profiles')
        qmissing = count_missing_data(qdf)
        rmissing = count_missing_data(rdf)

        total_samples = len(qdf) + len(rdf)
        missing_threshold = int(missing_threshold * total_samples)

        #Identify cols to remove

        for col in qmissing:
            count = qmissing[col]
            if not col in rmissing:
                rmissing[col] = 0
            rmissing[col] += count
            if rmissing[col] > missing_threshold:
                cols_to_remove.append(col)

        run_data['loci_removed'] = sorted(list(set(run_data['loci_removed']) | set(cols_to_remove)))

        if len(cols_to_remove) > 0:
            qdf = filter_columns(qdf, cols_to_remove)
            rdf = filter_columns(rdf, cols_to_remove)
    sys.stdout.flush()
    #convert profiles for fast dist calculations
    qlabels,qprofiles = convert_profiles(qdf)
    rlabels,rprofiles = convert_profiles(rdf)

    run_data['query_profile_info']['num_samples'] = len(qlabels)
    run_data['query_profile_info']['num_samples_pass'] = run_data['query_profile_info']['num_samples']
    run_data['ref_profile_info']['num_samples'] = len(qlabels)
    run_data['ref_profile_info']['num_samples_pass'] = run_data['ref_profile_info']['num_samples']

    if not skip:
        print(f'Performing QA/QC on input sample profiles')
        # Remove poor quality samples from the comparisons
        query_missing_data_counts = get_missing_loci_counts(qprofiles, qlabels, len(common_cols))
        ref_missing_data_counts = get_missing_loci_counts(rprofiles, rlabels, len(common_cols))
        query_samples_to_remove = flag_samples(query_missing_data_counts, sample_qual_thresh)
        run_data['query_profile_info']['failed_samples'] = query_samples_to_remove
        run_data['query_profile_info']['num_samples_pass'] = run_data['query_profile_info']['num_samples'] - len(query_samples_to_remove)
        ref_samples_to_remove = flag_samples(ref_missing_data_counts, sample_qual_thresh)
        run_data['ref_profile_info']['failed_samples'] = ref_samples_to_remove
        run_data['ref_profile_info']['num_samples_pass'] = run_data['ref_profile_info']['num_samples'] - len(ref_samples_to_remove)


        qlabels,qprofiles = filter_samples(qlabels, qprofiles, set(query_samples_to_remove) )
        rlabels, rprofiles = filter_samples(rlabels, rprofiles, set(ref_samples_to_remove))

        samples_to_remove = list(set(query_samples_to_remove) | set(ref_samples_to_remove))
        print(f'Prefilter Query size {len(qdf)}')
        qdf = qdf.drop(list(query_samples_to_remove))
        print(f'Postfilter Query size {len(qdf)}')
        rdf = rdf.drop(list(query_samples_to_remove))

    sys.stdout.flush()

    # write updated profiles
    print(f'Writting updated profiles to disk')
    write_profiles(qdf, os.path.join(outdir, f'query_profile.{file_type}'), file_type)
    run_data['query_profile_info']['parsed_file_path'] = os.path.join(outdir, f'query_profile.{file_type}')
    write_profiles(rdf, os.path.join(outdir, f'ref_profile.{file_type}'), file_type)
    run_data['ref_profile_info']['parsed_file_path'] = os.path.join(outdir, f'ref_profile.{file_type}')
    sys.stdout.flush()

    #Data frames no longer needed
    del(qdf)
    del (rdf)

    #Automatically determine batch size that fits in available memory
    num_records = len(qlabels) + len(rlabels)
    if len(qprofiles) == 0 or len(rprofiles) == 0:
        print(f'Error filtering parameters are too stringent there are 0 samples remaining: query:{len(qprofiles)} ref:{len(rprofiles)}')
        sys.exit()


    num_columns = len(qprofiles[0])
    byte_value_size = 8  #8 bytes for float64 which is the worst case
    if batch_size is None:
        batch_size = calc_batch_size(num_records,num_columns,byte_value_size,max_mem)
    print(f'Using a batch size of {batch_size}')
    sys.stdout.flush()

    #compute distances
    dist_matrix_file = os.path.join(outdir,f'dists.parquet')
    if os.path.isfile(dist_matrix_file):
        os.remove(dist_matrix_file)

    print(f'Calculating distances')
    if count_missing_sites:
        if dist_method == 'scaled':
            calc_distances_scaled_missing(qprofiles,qlabels,rprofiles,rlabels,dist_matrix_file,batch_size)
        else:
            calc_distances_hamming_missing(qprofiles, qlabels, rprofiles, rlabels, dist_matrix_file,batch_size)
    else:
        if dist_method == 'scaled':
            calc_distances_scaled(qprofiles,qlabels,rprofiles,rlabels,dist_matrix_file,batch_size)
        else:
            calc_distances_hamming(qprofiles, qlabels, rprofiles, rlabels, dist_matrix_file,batch_size)

    sys.stdout.flush()

    #format output for output format
    results_file = os.path.join(outdir,"results.{}".format(file_type))
    run_data['result_file'] = results_file
    write_dist_results(dist_matrix_file,
                       results_file, outfmt,
                       file_type, batch_size=batch_size, threshold=match_threshold)

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    sys.stdout.flush()

    with open(os.path.join(outdir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))

    if os.path.isfile(dist_matrix_file):
        os.remove(dist_matrix_file)

    sys.stdout.flush()

# call main function
if __name__ == '__main__':
    main()
