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
   rounds = '4'                 #Sets number of racon rounds
   cutoff = '1250'              #Sets filtlong read length cutoff
   keep = '90'                  #Sets filtlong read percentage to keep
   qual = '12'                  #Sets filtlong quality weighting
   pairing = '-p'               #Sets BWA flag to -p using interleaved illumina reads. If IlumnaB used then BWA -p removed.
   n = int(1)                   #Used for racon round loop for output files.
   dirs = '/'
   bugfix = 'export LC_ALL=C; unset LANGUAGE' #fix for racon compiler issue on cluster. 


	
   try:
      opts, args = getopt.getopt(argv,"hi:r:l:m:g:t:o:u:c:k:q:p:",["idir=", "reads=", "ilumnA=", "ilumnB=", "genomesize=", "threds=", "odir=", "rounds=", "cutoff=","keep=","qual=","pairing="])
   except getopt.GetoptError:
      print ('wrench.py -r [path/to/ONT.fastq] -l [path to ILLUMINA1.fastq] -g [genome size estimate] -o [output prefex] \n'
             '\n'
             '-t <threads>                (default 6) \n'
             '-u <racon rounds>           (default 4) \n'
             '-c <read length cutoff>     (default 1250) \n'
             '-q <read quality weighting> (default 12) \n'
             '-k <percent reads to keep>  (default 90) ')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('wrench.py \n' 
                '-r --reads            <path/to/ONT.fastq> \n'
                '-l --ilumnA           <path to ILLUMINA_1.fastq> \n'
                '-g --genomesize       <genome size estimate> \n'
                '-o --odir             <output tag> \n'
                '-m --ilumnB           <optional: path to ILLUMNA_2.fastq> \n'
                '-t --threads    (int) <optional: number of threads> \n'
                '-u --rounds     (int) <optional: number of racon polishing rounds (default 4)> \n'
                '-c --cutoff     (int) <optional: read length < (int) discarded (default 1250) > \n'
                '-k --keep       (int) <optional: percent of reads to retain (default 90)> \n'
                '-q --qual       (int) <optional: quality weighting for ONT.fastq filtering (default 12) \n' 
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
                'Requirements and pipeline: \n'
                'Filtlong \n'
                'Flye \n '
                'Minimap2 \n'
                'Racon \n'
                'Medaka \n '
                'BWA \n'
                'Samtools \n'
                'Pilon \n ')

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

#Make dir
   os.mkdir(outputdir)
   os.mkdir(outputdir+dirs+polishing+outputdir)


#Read filtering
   filtlong = 'filtlong ' \
              '--min_length {0} ' \
              ' --keep_percent {1} ' \
              ' --mean_q_weight {2} ' \
              '{3} | gzip > {4}/reads.q.fastq.gz' .format(cutoff,
                                                          keep,
                                                          qual,
                                                          input_fastQ_ONT,
						          outputdir)
   os.system(filtlong)
   os.system(bugfix)
	
#Read assembly
   flye = 'flye ' \
          '--nano-raw {0}/reads.q.fastq.gz ' \
          '--plasmids ' \
          '--out-dir {0}/flye_assembly_{0} ' \
          '--genome-size {1} ' \
          '--threads {2}' .format(outputdir,
				  genomesize,
				  threads)
   os.system(flye)

#Racon polishing

   map1 = 'minimap2 ' \
          '-ax map-ont {0}/flye_assembly_{0}/assembly.fasta ' \
          '{0}/reads.q.fastq.gz > {0}/{1}/racon_mapped.sam' .format(outputdir,
								    polishing+outputdir)

   pol1 = 'racon ' \
          '-m 8 ' \
          '-x -6 ' \
          '-g -8 ' \
          '-w 500 ' \
          '-t {0} ' \
          '{1}/reads.q.fastq.gz ' \
          '{1}/{2}/racon_mapped.sam ' \
          '{1}/flye_assembly_{1}/assembly.fasta > {1}/{2}/racon_{3}.fasta' .format(threads,
										   outputdir,
										   polishing+outputdir,
										   n)

   os.system(map1)
   os.system(pol1)

#Additional racon rounds:

   for x in range (1,rounds):
      map_loop = 'minimap2 ' \
                 '-ax map-ont ' \
                 '{0}/{1}/racon_{2}.fasta ' \
                 'reads.q.fastq.gz > {0}/{1}/racon_mapped.sam' .format(outputdir,
								       polishing+outputdir,
								       n)
      pol_loop = 'racon ' \
                 '-m 8 ' \
                 '-x -6 ' \
                 '-g -8 ' \
                 '-w 500 ' \
                 '-t {0} ' \
                 'reads.q.fastq.gz ' \
                 '{1}/{2}/racon_mapped.sam ' \
                 '{1}/{2}/racon_{3}.fasta > {1}/{2}/racon_{4}.fasta' .format(threads,
									     outputdir,
									     polishing+outputdir,
									     n,
									     n + 1) #+1 to racon file name

      clean = 'rm ' \
              '{0}/{1}/racon_{2}.fasta' .format(outputdir,
						polishing+outputdir,
						n)
      os.system(map_loop)
      os.system(pol_loop)
      os.system(clean) #Removes previous racon round
      n = n + 1

#Medaka polishing
   med = 'medaka_consensus ' \
         '-i {0}/reads.q.fastq.gz ' \
         '-d {0}/{1}/racon_{2}.fasta ' \
         '-o {0}/medaka_{0} ' \
         '-t {3} ' \
         '-m r941_min_high_g303' .format(outputdir,
					 polishing+outputdir,
					 rounds,
					 threads)
   os.system(med)

#Pilon hybrid assembly

   bwai = 'bwa index ' \
          '{0}/medaka_{0}/consensus.fasta' .format(outputdir)

   bwaa = 'bwa mem ' \
          '-M ' \
          '{0} ' \
          '{1}/medaka_{1}/consensus.fasta ' \
          '{2} {3} ' \
          '| samtools view -hu -F 4 - ' \
          '| samtools sort - > {1}/medaka_{1}/{1}.out.bam' .format(pairing,
								   outputdir,
								   input_fastQ_IluminaA,
								   input_fastQ_IluminaB)

   sami = 'samtools index ' \
          '{0}/medaka_{0}/{0}.out.bam' .format(outputdir)

   pilon = 'pilon ' \
           '--genome {0}/medaka_{0}/consensus.fasta ' \
           '--frags {0}/medaka_{0}/{0}.out.bam ' \
           '--outdir results_{0}' .format(outputdir)
   
   os.system(bwai)
   os.system(bwaa)
   os.system(sami)
   os.system(pilon)

   #for opt in opts:
      #if opt == '-E':
          #print('Cleaning up files')
          #cleanup = 'rm -r {0}/medaka* {0}/racon* {0}/flye*' .format(outputdir)
          #os.system(cleanup)
          #print('Assembly complete')
          #sys.exit()
      #else
          #print('Assembly complete')
          #sys.exit()
   print('done')

if __name__ == "__main__":
   main(sys.argv[1:])
