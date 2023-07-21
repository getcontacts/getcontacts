#!/usr/bin/env python
"""
Computes a clustered time-segmentation from a multi-frame contact file using TICC [1].

The input is a contact-file and a desired number of clusters, k. The output is a file
with the contents of the segmentation. It can either be a tab-separated file where each
line records a frame-number and the assigned cluster.

    0   0
    1   0
    2   1
    3   1
    4   1
    ...

It can also be a multi-flare file for use in displaying differences on the flareplot
webpage (gpcrviz.github.io/flareplot).
"""

import contact_calc.argparsers as ap
import argparse
from contact_calc.transformations import *
import numpy as np
import logging
import ticc
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
from contact_calc.transformations import res_contacts, gen_counts, gen_frequencies


def run_ticc(input_data, cluster_number=range(2, 11), process_pool_size=10, window_size=1,
             lambda_param=[1e-2], beta=[0.01, 0.1, 0.5, 10, 50, 100, 500], max_iters=1000, threshold=2e-5,
             BIC_Iters=15, logging_level=logging.CRITICAL):
    """
    Required Parameters:
    -- input_data: see input_format below
    -- output_filename: the output file name to write the cluster assignments
    Optional Parameters: BIC
    For each of these parameters, one can choose to specify:
        - a single number: this value will be used as the parameter
        - a list of numbers: the solver will use grid search on the BIC to choose the parameter
        - not specified: the solver will grid search on a default range (listed) to choose the parameter
    -- cluster_number: The number of clusters to classify. Default: BIC on [2...10]
    -- lambda_param: sparsity penalty. Default: BIC on 11e-2]
    -- beta: the switching penalty. If not specified, BIC on [50, 100, 200, 400]
    Other Optional Parameters:
    -- input_dimensions: if specified, will truncated SVD the matrix to the given number of features
       if the input is a graph, or PCA it if it's a matrix
    -- BIC_iters: if specified, will only run BIC tuning for the given number of iterations
    -- process_pool_size: the number of processes to spin off for optimization. Default 1
    -- window_size: The size of the window for each cluster. Default 1
    -- maxIters: the maximum number of iterations to allow TICC to run. Default 1000
    -- threshold: the convergence threshold. Default 2e-5
    -- covariance_filename: if not None, write the covariance into this file
    -- file_type is the type of data file. the data file must
       be a comma separated CSV. the options are:
       -- "matrix": a numpy matrix where each column is a feature and each
          row is a time step
       -- "graph": an adjacency list with each row having the form:
          <start label>, <end label>, value
    -- delimiter is the data file delimiter
    """
    logging.basicConfig(level=logging_level)
    # if input_format == 'graph':
    #     input_data = retrieveInputGraphData(
    #         input_filename, input_dimensions, delim=delimiter)
    # elif input_format == "matrix":
    #     input_data = np.loadtxt(input_filename, delimiter=delimiter)
    #     if input_dimensions is not None and input_dimensions < np.shape(input_data)[1]:
    #         pca = PCA(n_components=input_dimensions)
    #         input_data = pca.fit_transform(input_data)

    # print("Data shape %s, %s" % (np.shape(input_data)[0], np.shape(input_data)[1]))

    # get params via BIC
    cluster_number = cluster_number if isinstance(cluster_number, list) else list(cluster_number)
    beta = beta if isinstance(beta, list) else [beta]
    lambda_param = lambda_param if isinstance(lambda_param, list) else [lambda_param]
    BIC_Iters = max_iters if BIC_Iters is None else BIC_Iters
    problem_instance = ticc.ProblemInstance(input_data=input_data, window_size=window_size,
                                            maxIters=BIC_Iters, threshold=threshold)
    clusterResults = ticc.runHyperParameterTuning(beta, lambda_param, cluster_number, process_pool_size,
                                                  problem_instance)
    final_results = []
    for cluster_number, resultPackage in clusterResults:
        params, results, score = resultPackage
        beta, lambda_param = params
        # print("Via BIC with score %s, using params beta: %s, clusterNum %s, lambda %s" % (
        #     score, beta, cluster_number, lambda_param))
        # perform real run
        if BIC_Iters == max_iters:  # already performed the full run
            (cluster_assignments, cluster_MRFs) = results
        else:
            (cluster_assignment, cluster_MRFs) = ticc.solve(
                window_size=window_size, number_of_clusters=cluster_number, lambda_parameter=lambda_param,
                beta=beta, maxIters=max_iters, threshold=threshold,
                input_data=input_data, num_processes=process_pool_size, logging_level=logging_level)
        # outstream = "%s_%s" % (cluster_number, output_filename)
        # np.savetxt(outstream, cluster_assignment, fmt='%d', delimiter=',')
        final_results.append(
            (cluster_assignment, cluster_MRFs, (beta, lambda_param, cluster_number)))
    return final_results


def featurize_contacts(residue_contacts, dimensions):
    mapping = {}  # edge to value
    sparse_cols = []  # list of indices that should be 1
    counter = 0
    curr_timestamp = None

    for contact in residue_contacts:
        timestamp = contact[0]
        key = "%s_%s" % (contact[1], contact[2])
        if timestamp != curr_timestamp:  # new time
            curr_timestamp = timestamp
            sparse_cols.append(set())
        if key not in mapping:  # a new feature
            # assign this key to the current counter value
            mapping[key] = counter
            counter += 1
        # assign this feature into the current time step
        sparse_cols[-1].add(mapping[key])

    num_cols = len(mapping.keys())
    if dimensions is None or num_cols <= dimensions:  # do not need to SVD
        rows = []
        for indices in sparse_cols:
            # indices is a set
            row = [1.0 if i in indices else 0.0 for i in range(num_cols)]
            rows.append(row)
        return np.array(rows)
    else:
        # need truncated SVD
        data = []
        rows = []
        cols = []
        for i, indices in enumerate(sparse_cols):  # row
            for j in range(num_cols):  # col
                if j in indices:
                    data.append(1)
                    rows.append(i)
                    cols.append(j)
        mat = csr_matrix((data, (rows, cols)), shape=(len(sparse_cols), num_cols))
        solver = TruncatedSVD(n_components=dimensions)
        return solver.fit_transform(mat)


def main(argv=None):
    # Parse arguments
    parser = ap.PrintUsageParser(__doc__)
    parser.add_argument("--input_contacts",
                        type=argparse.FileType('r'),
                        required=True,
                        metavar="FILE",
                        help="Path to contact file")
    parser.add_argument("--clusters",
                        type=int,
                        required=False,
                        nargs="+",
                        default=[2, 5, 10],
                        metavar="INT",
                        help="Number of clusters [default: 2 5 10]")
    parser.add_argument("--tab_output",
                        type=str,
                        required=False,
                        metavar="FILE",
                        help="Path to TICC output file (tab-separated time/cluster indicators)")
    parser.add_argument("--frequency_output",
                        type=str,
                        required=False,
                        metavar="FILE",
                        help="Prefix to TICC output files (one res-frequency file for each cluster)")
    parser.add_argument("--beta",
                        type=int,
                        required=False,
                        nargs="+",
                        default=[10, 50, 100],
                        metavar="INT",
                        help="Beta parameter [default: 10 50 100]")
    parser.add_argument("--max_dimension",
                        type=int,
                        required=False,
                        default=50,
                        metavar="INT",
                        help="Max number of dimensions [default: 50]")
    args = parser.parse_args(argv)

    # Check output format and call corresponding function(s)
    if all(a is None for a in [args.tab_output, args.frequency_output]):
        parser.error("--tab_output or --frequency_output must be specified")

    print("Reading atomic contacts from " + args.input_contacts.name)
    atomic_contacts, num_frames = parse_contacts(args.input_contacts)
    args.input_contacts.close()

    print("Converting atomic contacts to residue contacts")
    residue_contacts = res_contacts(atomic_contacts)

    print("Performing dimensionality reduction")
    time_matrix = featurize_contacts(residue_contacts, args.max_dimension)

    print("Running TICC (clustered time-segmentation)")
    segmentation = run_ticc(time_matrix, cluster_number=args.clusters, beta=args.beta)

    if args.tab_output is not None:
        print("Writing time-segments to " + args.tab_output)
        with open(args.tab_output, "w") as f:
            f.writelines(map(lambda l: str(int(l)) + "\n", segmentation[0][0]))

    if args.frequency_output is not None:
        k = segmentation[0][2][2]
        for c in range(k):
            cluster_frames = set([frame for frame, cluster in enumerate(segmentation[0][0]) if cluster == c])
            cluster_contacts = [contact for contact in residue_contacts if contact[0] in cluster_frames]
            num_frames = len(cluster_frames)

            counts = gen_counts(cluster_contacts)
            total_frames, frequencies = gen_frequencies([(num_frames, counts)])

            fname = "%s_resfreq_cluster%03d.tsv" % (args.frequency_output, c)
            print("Writing frequency-flare to " + fname)
            with open(fname, "w") as output_file:
                output_file.write('#\ttotal_frames:%d\tinteraction_types:all\n' % total_frames)
                output_file.write('#\tColumns:\tresidue_1,\tresidue_2\tframe_count\tcontact_frequency\n')
                for (res1, res2), (count, frequency) in frequencies.items():
                    output_file.write('\t'.join([res1, res2, "%.3f" % frequency]) + "\n")



if __name__ == "__main__":
    main()

__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
