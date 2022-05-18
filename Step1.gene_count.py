#####1.Run gene count analysis####3
#!/usr/bin/env python
import argparse
import os
import glob
import subprocess
import itertools
cpu = 3
####set the dir of hisat2, samtool and stringTie
hisat2 = '/software/hisat2-2.0.5/hisat2'
samtools ='/software/samtools-1.3.1/samtools'
stringTie = '/software/stringtie-1.2.4.Linux_x86_64/stringtie'
parser = argparse.ArgumentParser(description='')
parser.add_argument('-d', '--dataDir',
                       dest='dataDir',
                       action='store',
                       help='Directory containg  reads',
                       )
parser.add_argument('-s', '--Spe',
                       dest='Spe',
                       action='store',
                       help='Species of sample data',
                       )
parser.add_argument('-o', '--output',
                       dest='output',
                       action='store',
                       help='Ouput file dir',
                       )
args = parser.parse_args()
index={
    "Mouse":{
        "Reference":"mouse",#hisat2 index
        "gtf":"Mus_musculus.GRCm38.76.gtf"#annotation gtf
    },
    "Human": {
    "Reference": "Homo_sapiens_GRCh38",#hisat2 index
    "gtf": "human_gencode_vch38.gtf"#annotation gtf
}
}

def sort_bam(bamFile):
    key_prefix = os.path.splitext(os.path.split(bamFile)[1])[0].split('_')[0]
    SamtoolsCmd = '{0} sort {1} -o ./{2}_sorted.bam\n'.format(samtools, bamFile, key_prefix)
    return SamtoolsCmd

def expression(sortedBamFile,anno):
    key_prefix = os.path.splitext(os.path.split(sortedBamFile)[1])[0].split('_')[0]
    StringTieCmd = 'mkdir ./{2}\ncd ./{2}\n{0} ../{3} -p 1 -G {1} -o ./{2}.gtf -l {2} -A ./gene_abund.{2} -B -e\n'.format(stringTie, anno, key_prefix,
                                                                                          sortedBamFile)
    return StringTieCmd

readDir = os.path.abspath(args.dataDir)
spe=args.Spe

# get read pairs
reads = glob.glob(readDir + '*.gz')
visited, pairs = set(), {}
for fn in reads:
    prefix=os.path.split(fn)[-1].split('_')[0]
    if os.path.exists(readDir +prefix+'_1.fastq.gz') and os.path.exists(readDir +prefix+'_2.fastq.gz'):
        if prefix not in pairs:
            pairs[prefix]=[readDir +prefix+'_1.fastq.gz',readDir +prefix+'_2.fastq.gz']
single={}
for each in reads:
    if '_' in os.path.split(each)[-1]:
        prefix = os.path.split(each)[-1].split('_')[0]
    else:
        prefix = os.path.split(each)[-1].split('.')[0]
    if prefix not in pairs:
        single[prefix] = each

print('Find {} read pairs'.format(len(pairs)))
print('Find {} read single'.format(len(single)))

job_pair_num=0
for key, value in pairs.items():
    f1,f2=value
    job_pair_num += 1
    refPrefix = index[spe]["Reference"]
    anno = index[spe]["gtf"]
    for i in range(1, 9):
        indexFile = '{0}.{1}.ht2'.format(refPrefix, i)
        if not os.path.exists(indexFile) or not os.path.getsize(indexFile):
            print('cDNA HISAT2 index file: {0} does not exist or is zero in size'.format(indexFile))
            exit()
    cmd = '{7} -p {3} -x {0} -1 {1} -2 {2}| {4} view -bS - > {6}/{5}.bam\n'.format(
        refPrefix, f1, f2, cpu, samtools, key,readDir,hisat2)
    cur = os.getcwd()
    text1 = '#!/bin/bash\n#PBS -N {0}\n#PBS -l nodes=1:ppn=3\n#PBS -l walltime=1000:00:00\n#PBS -j oe\n#PBS -q silver\ncd {1}\n'.format(
        job_pair_num, cur)
    cmd1 = sort_bam(key + '.bam')
    cmd2 = expression(key + '_sorted.bam', anno)
    with open('{0}_exp.sh'.format(key), 'w') as af:
        af.write(text1)
        af.write(cmd)
        af.write(cmd1)
        af.write(cmd2)
job_single_num=0
for key, value in single.items():
    f1=value
    job_single_num += 1
    refPrefix = index[spe]["Reference"]
    anno = index[spe]["gtf"]
    for i in range(1, 9):
        indexFile = '{0}.{1}.ht2'.format(refPrefix, i)
        if not os.path.exists(indexFile) or not os.path.getsize(indexFile):
            print('cDNA HISAT2 index file: {0} does not exist or is zero in size'.format(indexFile))
            exit()
    cmd = '{6} -p {2} -x {0} -U {1} | {3} view -bS - > {5}/{4}.bam\n'.format(
        refPrefix, f1, cpu, samtools, key,readDir,hisat2)
    cur = os.getcwd()
    text1 = '#!/bin/bash\n#PBS -N {0}\n#PBS -l nodes=1:ppn=3\n#PBS -l walltime=1000:00:00\n#PBS -j oe\n#PBS -q silver\ncd {1}\n'.format(
        job_single_num, cur)
    cmd1 = sort_bam(key + '.bam')
    cmd2 = expression(key + '_sorted.bam', anno)
    with open('{0}_exp.sh'.format(key), 'w') as af:
        af.write(text1)
        af.write(cmd)
        af.write(cmd1)
        af.write(cmd2)