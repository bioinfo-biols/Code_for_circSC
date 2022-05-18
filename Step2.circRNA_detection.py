#####2.Run CIRI2+CIRI-AS analysis####
#!/usr/bin/env python
import argparse
import os
import glob
import subprocess
import itertools
ciri = '/software/CIRI_v2.0.6/CIRI2.pl'
bwa='/software/bwa-0.7.12/bwa'
cirias = '/software/CIRI_AS_v1.2.pl'
parser = argparse.ArgumentParser(description='Performing CIRI on the paired-end reads on the given directory')
parser.add_argument('-d', '--dataDir',
                    dest='dataDir',
                    action='store',
                    help='Directory containg paired-end reads')
parser.add_argument('-q', '--quene',
                    dest='quene',
                    action='store',
                    help='quene')
args = parser.parse_args()

index={
    "Homo sapiens":{
        "rRNAReference":"/database/hisat2_index/human_rRNA",
        "Reference":/database/hisat2_index/Homo_sapiens_GRCh38",
        "bwaindex":"/database/bwa_index/Homo_sapiens_GRCh38.fa",
        "gtf":"/annotation/human_gencode_vch38.gtf"
    },
    "Mus musculus":{
        "rRNAReference":"/database/hisat2_index/mouse_rRNA",
        "Reference":"/database/hisat2_index/mouse",
        "bwaindex":"/database/bwa_index/Mus_musculus.GRCm38.dna.primary_assembly.fa",
        "gtf":"/annotation/Mus_musculus.GRCm38.76.gtf"
    }
}
human_sample=[]
mouse_sample=[]
with open('all_download_data','r') as an:
    for line in an:
        aa = line.rstrip().split('\t')
        if aa[1] =="Homo sapiens":
            human_sample.append(aa[0])
        else:
            mouse_sample.append(aa[0])

readDir = os.path.abspath(args.dataDir)
# get read pairs
reads = glob.glob(readDir + '/*_*.fastq.gz')
reads1 = glob.glob(readDir + '/*.fastq.gz')
visited, pairs = set(), {}
for fn in reads:
    prefix=os.path.split(fn)[1].split('_')[0]
    if prefix not in pairs:
        pairs[prefix]=[readDir +'/'+prefix+'_1.fastq.gz',readDir +'/'+prefix+'_2.fastq.gz']
single={}
for each in reads1:
    prefix = os.path.split(each)[1].split('.')[0].split('_')[0]
    if prefix not in pairs:
        single[prefix] = each
current_dir = os.getcwd()
job_pair_num = 0
for key, value in pairs.iteritems():
    if key in human_sample:
        refPrefix = index["Homo sapiens"]["bwaindex"]
        gtf = os.path.abspath(index["Homo sapiens"]["gtf"])
    else:
        refPrefix = index["Mus musculus"]["bwaindex"]
        gtf = os.path.abspath(index["Mus musculus"]["gtf"])
    for i in ['amb', 'ann', 'bwt', 'pac', 'sa']:
        indexFile = '{0}.{1}'.format(refPrefix, i)
        if not os.path.exists(indexFile) or not os.path.getsize(indexFile):
            print('cDNA BWA index file: {0} does not exist or is zero in size'.format(indexFile))
            exit()
    job_pair_num += 1
    f1, f2 = value
    cmd1 = '#!/bin/bash\n#PBS -N {0}\n#PBS -l nodes=1:ppn=6\n#PBS -l walltime=1000:00:00\n#PBS -j oe\n#PBS -q {2}\ncd {1}\n'.format(job_pair_num, current_dir, args.quene)
    cmd2 = 'bwa mem -T 19 -t 16 {0} {1} {2}> {3}.sam\n'.format(refPrefix, f1,f2, key)
    cmd3 = 'perl {0} -I {1}.sam -O {1}.ciri -F {2} -A {3} -T 16 -0'.format(ciri, key, refPrefix, gtf)
    cmd4 = 'perl {0} -S {1}.sam -C {1}.ciri -F {2} -A {3} -O {1} -D yes'.format(cirias, k, refPrefix,
                                                                                              gtf)
    with open('{0}.sh'.format(key), 'w') as af:
        af.write(cmd1)
        af.write(cmd2)
        af.write(cmd3)
        af.write(cmd4)

job_single_num = 0
for key, value in single.iteritems():
    job_single_num += 1
    if key in human_sample:
        refPrefix = index["Homo sapiens"]["bwaindex"]
        gtf = os.path.abspath(index["Homo sapiens"]["gtf"])
    else:
        refPrefix = index["Mus musculus"]["bwaindex"]
        gtf = os.path.abspath(index["Mus musculus"]["gtf"])
    for i in ['amb', 'ann', 'bwt', 'pac', 'sa']:
        indexFile = '{0}.{1}'.format(refPrefix, i)
        if not os.path.exists(indexFile) or not os.path.getsize(indexFile):
            print('cDNA BWA index file: {0} does not exist or is zero in size'.format(indexFile))
            exit()
    f1 = value
    cmd1 = '#!/bin/bash\n#PBS -N {0}\n#PBS -l nodes=1:ppn=6\n#PBS -l walltime=1000:00:00\n#PBS -j oe\n#PBS -q {2}\ncd {1}\n'.format(job_single_num, current_dir, args.quene)
    cmd2 = 'bwa mem -T 19 -t 16 {0} {1} > {2}.sam\n'.format(refPrefix, f1, key)
    cmd3 = 'perl {0} -I {1}.sam -O {1}.ciri -F {2} -A {3} -T 16 -0'.format(ciri, key, refPrefix, gtf)
    with open('{0}.sh'.format(key), 'w') as af:
        af.write(cmd1)
        af.write(cmd2)
        af.write(cmd3)