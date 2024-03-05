import warnings
import numpy as np
from collections import namedtuple
from scipy.signal import find_peaks_cwt
from matplotlib import pyplot as plt
import sys
import os
import argparse

import matplotlib as mpl
mpl.use('Agg')

warnings.filterwarnings("ignore")

QCResult = namedtuple('QCResult', ['metric', 'qc_pass', 'message'])
INF = float("inf")

class QCCheck(object):
    def __init__(self, metric):
        self.metric = metric

    def check(self, value):
        return True

    def message(self, value, qc_pass):
        return ('{}\tOK'.format(value) if qc_pass
                else '{}\tFailed'.format(value))

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, qc_pass, self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        return ('{}\tOK'.format(value) if qc_pass else
                '{}\tout of range [{}, {}]'.format(value, self.lower,
                                                  self.upper))


class QCLessThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, upper):
        super(QCLessThanEqualCheck, self).__init__(metric, -INF, upper)


class QCGreaterThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, lower):
        super(QCGreaterThanEqualCheck, self).__init__(metric, lower, INF)


class QCHasElementInRange(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCHasElementInRange, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, elems):
        return (len([elem for elem in elems
                    if self.lower <= elem <= self.upper]) > 0)

    def message(self, elems, qc_pass):
        return ('OK' if qc_pass else
                'Cannot find element in range [{}, {}]'.format(
                    self.lower, self.upper))


def read_picard_histogram(data_file):
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## HISTOGRAM'):
                break
        data = np.loadtxt(fp, skiprows=1)

    return data

def get_insert_distribution(final_bam, prefix):
    '''
    Calls Picard CollectInsertSizeMetrics
    '''
    logging.info('insert size distribution...')
    insert_data = '{0}.inserts.hist_data.log'.format(prefix)
    insert_plot = '{0}.inserts.hist_graph.pdf'.format(prefix)
    graph_insert_dist = ('java -Xmx6G -XX:ParallelGCThreads=1 -jar '
                         '{3} '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'VERBOSITY=ERROR QUIET=TRUE '
                         'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot,
                                                             locate_picard())
    logging.info(graph_insert_dist)
    os.system(graph_insert_dist)
    return insert_data, insert_plot

def fragment_length_qc(data, prefix):
    results = []

    NFR_UPPER_LIMIT = 150
    MONO_NUC_LOWER_LIMIT = 150
    MONO_NUC_UPPER_LIMIT = 300

    # % of NFR vs res
    nfr_reads = data[data[:,0] < NFR_UPPER_LIMIT][:,1]
    percent_nfr = nfr_reads.sum() / data[:,1].sum()
    results.append(
        QCGreaterThanEqualCheck('Fraction of reads in NFR', 0.4)(percent_nfr))

    # % of NFR vs mononucleosome
    mono_nuc_reads = data[
        (data[:,0] > MONO_NUC_LOWER_LIMIT) &
        (data[:,0] <= MONO_NUC_UPPER_LIMIT)][:,1]
    
    percent_nfr_vs_mono_nuc = (
        nfr_reads.sum() /
        mono_nuc_reads.sum())
    results.append(
        QCGreaterThanEqualCheck('NFR / mono-nuc reads', 2.5)(
            percent_nfr_vs_mono_nuc))

    # peak locations
    pos_start_val = data[0,0] # this may be greater than 0
    peaks = find_peaks_cwt(data[:, 1], np.array([25]))
    nuc_range_metrics = [('Presence of NFR peak', 20 - pos_start_val, 90 - pos_start_val),
                         ('Presence of Mono-Nuc peak', 120 - pos_start_val, 250 - pos_start_val),
                         ('Presence of Di-Nuc peak', 300 - pos_start_val, 500 - pos_start_val)]
    for range_metric in nuc_range_metrics:
        results.append(QCHasElementInRange(*range_metric)(peaks))

    out = prefix + '.nucleosomal.qc'
    with open(out, 'w') as fp:
        for elem in results:            
            fp.write('\t'.join([elem.metric, str(elem.qc_pass), elem.message]) + '\n')

    return out

def fragment_length_plot(data_file, prefix, peaks=None):
    try:
        data = read_picard_histogram(data_file)
    except IOError:
        return ''
    except TypeError:
        return ''

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

    if peaks:
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')

    # plot_img = BytesIO()
    # fig.savefig(plot_img, format='png')
    plot_png = prefix + '.fraglen_dist.png'
    fig.savefig(plot_png, format='png')
    
    print("{} saved!".format(plot_png))

    return plot_png


insert_data = sys.argv[1]
OUTPUT_PREFIX = sys.argv[2]

nucleosomal_qc = fragment_length_qc(read_picard_histogram(insert_data), OUTPUT_PREFIX)
fraglen_dist_plot = fragment_length_plot(insert_data, OUTPUT_PREFIX)




