from .fastq import Fastq
import numpy as np
import pandas as pd
import statistics
import os

EXE = "java -jar /usr/local/bin/picard.jar"

def cal_insert_picard(base_name, ref, adapters=None):
    print(base_name)
    input_sam = "{}.sam".format(base_name)
    output_sam = "{}_sorted.sam".format(base_name)
    report = "{}.report".format(base_name)
    hist = "{}.hist.pdf".format(base_name)
    output_alignment_summary = "{}.alignment.summary".format(base_name)
    output_gc_metrices = "{}.gc.metrices".format(base_name)
    output_gc_chart = "{}.gc.chart.pdf".format(base_name)
    output_gc_summary = "{}.gc.summary".format(base_name)
    os.system("{} SortSam I={} O={} SO=coordinate".format(EXE, input_sam, output_sam))
    os.system("{} CollectInsertSizeMetrics I={} O={} H={}".format(EXE, output_sam, report, hist))
    if adapters is not None:
        os.system("{} CollectAlignmentSummaryMetrics R={} I={} O={} ADAPTER_SEQUENCE={} ADAPTER_SEQUENCE={}".format(EXE, ref, output_sam, output_alignment_summary, adapters[0], adapters[1]))
    else:
        os.system("{} CollectAlignmentSummaryMetrics R={} I={} O={}".format(EXE, ref, output_sam, output_alignment_summary))
    os.system("{} CollectGcBiasMetrics I={} O={} CHART={} S={} R={}".format(EXE, output_sam, output_gc_metrices, output_gc_chart, output_gc_summary, ref))

    return output_sam

def cal_real_distribution(f_name):
    f = open(f_name, "r")
    for i in range(34):
        f.readline()
    header = ['# readId', 'position', 'insertSize']
    rows = []
    for i in range(10000):
        strs = f.readline()[:-1].split('\t')
        row = [strs[0], int(strs[3]), int(strs[5])]
        rows.append(row)
        
    df = pd.DataFrame(rows, columns=header)
    print(df['insertSize'].mean())
    print(df['insertSize'].std())
    df.to_excel("Sim_100_100.xlsx")

if __name__ == "__main__":
    pass