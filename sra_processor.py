import os
import subprocess

sra_list = os.listdir('sra')[1::2]

os.chdir('..')

sra_command = 'sra_toolkit/bin/./fastq-dump'

for sra in sra_list:
    sra_name = sra.split('.')[0]
    print 'Processing SRA file {}'.format(sra_name)
    subprocess.call([sra_command + ' marko_rna/sra/{}'.format(sra) + ' marko_rna/sra_fastq/{}.fastq'.format(sra_name)], shell=True)

print 'Done!'
