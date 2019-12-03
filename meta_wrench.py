#!/usr/bin/python

import sys
import os
import getopt


#Define current working dir
cwd = os.getcwd()

#Define inputs
def main(argv):
   input_fastQ_ONT = ''
   outputdir = ''
   genomesize = ''              #Estimate of genome size for assembly
   input_fastQ_IluminaA = ''    #Input for interleaved paired ends and primary illumina reads
   input_fastQ_IluminaB = ''    #Input for secondary illumina reads
   threads = ''                 #Sets threads to use for all sections
   polishing = 'racon_'         #prefex for polishing directory naming
   rounds = ''                  #Sets number of racon rounds
   cutoff = '1500'              #Sets filtlong read length cutoff
   keep = '90'                  #Sets filtlong read percentage to keep
   qual = '15'                  #Sets filtlong quality weighting
   pairing = '-p'               #Sets BWA flag to -p using interleaved illumina reads. If IlumnaB used then BWA -p removed.
   n = int(1)                   #Used for racon round loop for output files.

   try:
      opts, args = getopt.getopt(argv,"hi:r:l:m:g:t:o:u:c:k:q:p:",["idir=", "reads=", "ilumnA=", "ilumnB=", "genomesize=", "threds=", "odir=", "rounds=", "cutoff=","keep=","qual=","pairing="])
   except getopt.GetoptError:
      print ('wrench.py -r [path/to/ONT.fastq] -l [path to ILLUMINA1.fastq] -g [genome size estimate] -u [racon rounds] -o [output prefex] \n'
             '\n'
             '-t <threads>                (default 6) \n'
             '-u <racon rounds>           (default 4)\n'
             '-c <read length cutoff>     (default 1500bp)\n'
             '-q <read quality weighting> (default 15)\n'
             '-k <percent reads to keep>  (default 90) ')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('wrench.py \n' 
                '-r --reads            <full path/to/ONT.fastq> \n'
                '-l --ilumnA           <full path to ILLUMINA_1.fastq> \n'
                '-m --ilumnB           <optional: full path to ILLUMNA_2.fastq> \n'
                '-g --genomesize       <genome size estimate> \n'
                '-o --odir             <output tag> \n'
                '-t --threads    (int) <optional: number of threads> \n'
                '-u --rounds     (int) <number of racon polishing rounds: 4 recommended> \n'
                '-c --cutoff     (int) <read length: < (int) discarded > \n'
                '-k --keep       (int) <percent of reads to retain with q weighting> \n'
                '-q --qual       (int) <quality weighting for ONT.fastq filtering \n' 
                '\n'
                'wrench.py -r [path/to/ONT.fastq] -l [path to ILLUMINA1.fastq] -g [genome size estimate] -o [output prefex] \n'
                '\n'
                'OPTIONAL \n'
                '-m <path to illumina B reads>'
                '-t <threads>                (default 6) \n'
                '-u <racon rounds>           (default 4)\n'
                '-c <read length cutoff>     (default 1500bp)\n'
                '-q <read quality weighting> (default 15)\n'
                '-k <percent reads to keep>  (default 90) \n'
                '\n'
                'Requirements \n'
                'Filtlong \n'
                'Flye \n '
                'Minimap2 \n'
                'Racon \n'
                'Medaka \n '
                'Pilon \n '
                'BWA \n'
                'samtools')

         sys.exit()
      elif opt in ("-i", "--idir"):
         input_fast5_ONT = arg
      elif opt in ("-r", "--reads"):
         input_fastQ_ONT = arg
      elif opt in ("-l", "--ilumnA"):
         input_fastQ_IluminaA = arg
      elif opt in ("-m", "--ilumnB"):
         input_fastQ_IluminaB = arg
         pairing = ''
      elif opt in ("-o", "--odir"):
         outputdir = arg
      elif opt in ("-g", "--genomesize"):
         genomesize = arg
      elif opt in ("-t", "--threds"):
         threads = arg
      elif opt in ("-u", "--rounds"):
         rounds = int(arg)
      elif opt in ("-c", "--cutoff"):
         cutoff = int(arg)
      elif opt in ("-k", "--keep"):
         keep = int(arg)
      elif opt in ("-q", "--qual"):
         qual = int(arg)

#Make output directorys
   os.mkdir(outputdir)
   os.chdir(outputdir)
   os.mkdir(polishing+outputdir)


#Read filtering
   filtlong = 'filtlong ' \
              '--min_length {} ' \
              ' --keep_percent {} ' \
              ' --mean_q_weight {} ' \
              ' {} | gzip > reads.q.fastq.gz' .format(cutoff,
                                                      keep,
                                                      qual,
                                                      input_fastQ_ONT)
   os.system(filtlong)

#Read assembly
   flye = 'flye ' \
          '--nano-raw reads.q.fastq.gz ' \
          '--meta ' \
          '--out-dir flye_assembly_{} ' \
          '--genome-size {} ' \
          '--threads {}' .format(outputdir,
                                 genomesize,
                                 threads)
   os.system(flye)

#Racon polishing

   map1 = 'minimap2 ' \
          '-ax map-ont flye_assembly_{}/assembly.fasta ' \
          'reads.q.fastq.gz > {}/racon_mapped.sam' .format(outputdir,
                                                           polishing+outputdir)

   pol1 = 'racon ' \
          '-m 8 ' \
          '-x -6 ' \
          '-g -8 ' \
          '-w 500 ' \
          'reads.q.fastq.gz ' \
          '{}/racon_mapped.sam ' \
          'flye_assembly_{}/assembly.fasta > {}/racon_{}.fasta' .format(polishing+outputdir,
                                                                        outputdir,
                                                                        polishing+outputdir,
                                                                        n)

   os.system(map1)
   os.system(pol1)

#Additional racon rounds:

   for x in range (1,rounds):
      map_loop = 'minimap2 ' \
                 '-ax map-ont ' \
                 '{}/racon_{}.fasta ' \
                 'reads.q.fastq.gz > {}/racon_mapped.sam' .format(polishing+outputdir,
                                                                  n,
                                                                  polishing+outputdir)
      pol_loop = 'racon ' \
                 '-m 8 ' \
                 '-x -6 ' \
                 '-g -8 ' \
                 '-w 500 ' \
                 'reads.q.fastq.gz ' \
                 '{}/racon_mapped.sam ' \
                 '{}/racon_{}.fasta > {}/racon_{}.fasta' .format(polishing+outputdir,
                                                                 polishing+outputdir,
                                                                 n,
                                                                 polishing+outputdir,
                                                                 n + 1) #+1 to racon file name

      clean = 'rm ' \
              '{}/racon_{}.fasta' .format(polishing+outputdir,
                                          n)
      os.system(map_loop)
      os.system(pol_loop)
      os.system(clean) #Removes previous racon round
      n = n + 1

#Medaka polishing
   med = 'medaka_consensus ' \
         '-i reads.q.fastq.gz ' \
         '-d {}/racon_{}.fasta ' \
         '-o medaka_{} ' \
         '-t {} ' \
         '-m r941_min_high_g303' .format(polishing+outputdir,
                                         rounds,
                                         outputdir,
                                         threads)
   os.system(med)

#Pilon hybrid assembly

   bwai = 'bwa index ' \
          './medaka_{}/consensus.fasta' .format(outputdir)

   bwaa = 'bwa mem ' \
          '-M ' \
          '{} ' \
          'medaka_{}/consensus.fasta ' \
          '{} {} ' \
          '| samtools view -hu -F 4 - ' \
          '| samtools sort - > medaka_{}/{}.out.bam' .format(pairing,
                                                             outputdir,
                                                             input_fastQ_IluminaA,
                                                             input_fastQ_IluminaB,
                                                             outputdir,
                                                             outputdir)

   sami = 'samtools index ' \
          'medaka_{}/{}.out.bam' .format(outputdir,
                                         outputdir)

   pilon = 'pilon ' \
           '--genome medaka_{}/consensus.fasta ' \
           '--frags medaka_{}/{}.out.bam ' \
           '--outdir results_{}' .format(outputdir,
                                         outputdir,
                                         outputdir,
                                         outputdir)
   os.system(bwai)
   os.system(bwaa)
   os.system(sami)
   os.system(pilon)

   #for opt in opts:
      #if opt == '-E':
          #print('Cleaning up files')
          #cleanup = 'rm -r medaka* racon* flye*'
          #os.system(cleanup)
          #print('Assembly complete')
          #sys.exit()
      #else
          #print('Assembly complete')
          #sys.exit()
   print('done')

if __name__ == "__main__":
   main(sys.argv[1:])
