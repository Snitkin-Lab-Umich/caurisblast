

#from Bio import SeqIO
import os
import subprocess
import argparse

# def file_reader(input_file_path):
#     genelist = []
#     for record in SeqIO.parse(input_file_path):
#         genelist += ['>' + record.id,record.seq]

def make_query(dir_path,outputfile):
    clist = []
    flist = os.listdir(dir_path)
    for f in flist:
        if sum([x in f for x in ['.txt','.fasta','.fa']]) == 0 and f != outputfile:
            print(f'Skipping {f}')
            continue
        #print(f'Adding {f} to query file')
        clist += [f]
    with open(outputfile, 'w') as fhout:
        #subprocess.call(['cat'] + clist, 
        #                cwd = dir_path, stdout=output, stderr=debug)
        for fname in clist:
            with open(dir_path+fname,'r') as fhin:
                print(f'Adding {fname} to query file')
                for line in fhin:
                    _ = fhout.write(line)
                if line[-1] != '\n':
                    _ = fhout.write('\n')
                # this should ensure that each individual fasta ends with a new line character as it's added to the combined query
    return(outputfile)

# input_path must be a single path
# if a directory is provided, combine the files first
def check_query(input_path,output_name):
    if os.path.isdir(input_path):
        query_fasta = make_query(input_path + '/',output_name + '_combined_query.fasta')
    elif not os.path.isfile(input_path):
        print(f'Cannot locate file or directory at {input_path}')
        quit(1)
    else:
        query_fasta = input_path
    return(query_fasta)

def make_subject_fasta(input_dir,output_name):
    flist = [x for x in os.listdir(input_dir) if x.endswith('.fasta') or x.endswith('.fa')]
    if len(flist) < 1:
        print('Cannot find .fasta or .fa files for BLAST database in the provided directory')
        quit(1)
    outpath = 'db/' + output_name + '_blastdb.fasta'
    with open(outpath,'w') as fhout:
        for fname in flist:
            prefix = fname.split('.fa')[0]
            with open(input_dir + fname, 'r') as fhin:
                for line in fhin:
                    if line.startswith('>'):
                        line = '>' + prefix + '_' + line.split('>')[1]
                    _ = fhout.write(line)
                if line[-1] != '\n':
                    _ = fhout.write('\n')
    return(outpath)

def check_subject(input_path,output_name,debuglog = 'logs/debug_log.txt'):
    if os.path.isdir(input_path):
        subject_fasta = make_subject_fasta(input_path + '/',output_name)
    elif os.path.isfile(input_path) and (input_path.endswith('.fasta') or input_path.endswith('.fa')):
        subject_fasta = 'db/' + output_name + '_blastdb.fasta'
        with open(debuglog, 'a') as debug:
            command = ['cp',input_path,subject_fasta]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
    # with either method, the db/ directory should contain the .fasta file needed for making the database
    else:
        print('Cannot find .fasta file for BLAST database')
        quit(1)
    return(subject_fasta)

def make_blast_db(input_fasta,output_name,type = 'nucl',container = None,debuglog = 'logs/debug_log.txt'):
    with open(debuglog, 'a') as debug:
        if container is None:
            command = ['makeblastdb','-in',input_fasta,'-title',output_name,'-dbtype',type]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
        else:
            command = ['singularity','exec',container,'makeblastdb','-in',input_fasta,'-title',output_name,'-dbtype',type]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
        return(input_fasta)

def run_blast(input,database,output,type = 'nucl',format = '10',container = None,debuglog = 'logs/debug_log.txt'):
    if type == 'nucl':
        blastype = 'blastn'
    elif type == 'prot':
        blastype = 'blastp'
    else:
        print('BLAST type not recognized, please specify "nucl" or "prot"')
        quit(1)
    with open(debuglog, 'a') as debug:
        if container is None:
            command = [blastype,'-db',database,'-query',input,'-out',output,'-outfmt',format]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
        else:
            command = ['singularity','exec',container,blastype,'-db',database,'-query',input,'-out',output,'-outfmt',format]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
    print('Finished running BLAST!')
    print(f'Results file: {output}')
    print(f'Database location: {database}')
    print(f'Debug log and details: {debuglog}')

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--query','-q',type=str,
        help='''Provide a query sequence in fasta format (.fasta or .fa). This can be either a single file or a directory of files.''',
        default=None,required=True
        )
    parser.add_argument(
        '--subject','-s',type=str,
        help='''Provide a subject to perform the search on. This can be a single file or a directory of files, all in the fasta format. 
        A new BLAST database will be generated from the provided sequences.
        If you are using an existing database, provide the path to the database and use the "--database" flag.''',
        default=None,required=True
        )
    parser.add_argument(
        '--database','-d',action='store_true',help='''Add this flag if you provided a pre-built BLAST database with "--subject".'''
        )
    parser.add_argument(
        '--type','-t',type=str,choices=['n','nucl','nucleotide','p','prot','protein'],
        help='''Specify whether this should be a nucleotide BLAST or a protein BLAST.
        Note that the database must match the type of BLAST you perform.''',
        default='nucl',required=True
        )
    parser.add_argument(
        '--name','-n',type=str,help='''(Optional) Provide a name for the run. This name will be used for all database and search outputs.''',
        default='new_search_'
        )
    parser.add_argument(
        '--format','-f',type=str,help='''(Optional) Provide an alternate format for the BLAST search output. The default is format 10.''',
        default='10'
        )
    parser.add_argument(
        '--annotation','-a',type=str,help='''(Optional) Provide a protein annotation file for use with a nucleotide BLAST.
        An additional output will be created that indicates which proteins your query sequences appear in, if any. 
        This assumes the .fasta you used for your BLAST database has headers formatted as: ">[sample]_[sequence name]".''',
        default=None
        )
    parser.add_argument(
        '--container','-c',type=str,help='''(Optional) Provide a singularity container for the BLAST+ software.
        The specified container will be used when running the BLAST, instead of assuming you have BLAST+ installed locally.''',
        default=None
        )
    parser.add_argument(
        '--repeatmasker','-r',type=str,help='''(Optional) Provide a path to a RepeatMasker container file. This will be used to 
        mask repetitive sequences in the query and subject.''',
        default=None
        )
    args = parser.parse_args()
    # handle any missing or incorrect arguments
    if args.type in ['n','nucl','nucleotide']:
        args.type = 'nucl'
    elif args.type in ['p','prot','protein']:
        args.type = 'prot'
    else:
        print(f'Unrecognized BLAST search type "{args.type}"')
        quit(1)
    if args.query is None or args.subject is None:
        print('Missing query or subject sequence')
        quit(1)
    with open('logs/debug_log.txt','w') as fh:
        _ = fh.write(f'debug log for {args.name}\n')
    # change working directory to location of script
    for f in [args.query,args.subject,args.annotation]:
        if f is not None:
            if not os.path.exists(f):
                print(f'Could not locate file or directory at {f}')
                quit(1)
    args.query,args.subject,args.annotation = [os.path.abspath(x) if x is not None else x for x in [args.query,args.subject,args.annotation]]
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # generate the files and directories needed for the BLAST search in this directory
    query_fasta = check_query(args.query,args.name)
    if args.repeatmasker is not None:
        print('REPEAT MASKING CODE IN PROGRESS')
        # perform repeat masking on the query and subject
    if not args.database:
        subject_fasta = check_subject(args.subject,args.name)
        subject_fasta = make_blast_db(input_fasta = subject_fasta,output_name = args.name,type = args.type,container = args.container)
        # this will generate a blast database, if needed
    else:
        subject_fasta = args.subject
    search_output_name = 'results/' + args.name + '_' + args.type + '_blast_results.csv'
    run_blast(input=query_fasta, database=subject_fasta,output=search_output_name,type=args.type,format=args.format,container=args.container)
    
    

if __name__ == "__main__":
    main()







