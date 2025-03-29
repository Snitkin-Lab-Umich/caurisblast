

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
# either way, ensure the final query is located at temp/{output_name}/
def check_query(input_path,output_name,debuglog):
    final_query = f'results/{output_name}/{output_name}_query.fasta'
    if os.path.isdir(input_path):
        final_query = make_query(input_path + '/',final_query)
    elif not os.path.isfile(input_path):
        print(f'Cannot locate file or directory at {input_path}')
        quit(1)
    else:
        cdbg([['cp',input_path,final_query]],debuglog)
    return(final_query)

def make_subject_fasta(input_dir,outputfile):
    flist = [x for x in os.listdir(input_dir) if x.endswith('.fasta') or x.endswith('.fa')]
    if len(flist) < 1:
        print('Cannot find .fasta or .fa files for BLAST database in the provided directory')
        quit(1)
    #outpath = f'db/{output_name}' + output_name + '_blastdb.fasta'
    with open(outputfile,'w') as fhout:
        for fname in flist:
            prefix = fname.split('.fa')[0]
            with open(input_dir + fname, 'r') as fhin:
                for line in fhin:
                    if line.startswith('>'):
                        line = '>' + prefix + '_' + line.split('>')[1]
                    _ = fhout.write(line)
                if line[-1] != '\n':
                    _ = fhout.write('\n')
    return(outputfile)

def check_subject(input_path,output_name,debuglog):
    final_subject = f'db/{output_name}/{output_name}_blastdb.fasta'
    if os.path.isdir(input_path):
        final_subject = make_subject_fasta(input_path + '/',final_subject)
    elif os.path.isfile(input_path) and (input_path.endswith('.fasta') or input_path.endswith('.fa')):
        cdbg([['cp',input_path,final_subject]],debuglog)
    else:
        print('Cannot find .fasta file for BLAST database')
        quit(1)
    return(final_subject)

def make_blast_db(input_fasta,output_name,type = 'nucl',container = None,debuglog = 'logs/debug_log.txt',mask_data = None):
    with open(debuglog, 'a') as debug:
        if container is None:
            command = ['makeblastdb','-in',input_fasta,'-title',output_name,'-dbtype',type]
            if mask_data is not None:
                command+=['-mask_data',mask_data]
            cdbg([command],debuglog)
        else:
            command = ['singularity','exec',container,'makeblastdb','-in',input_fasta,'-title',output_name,'-dbtype',type]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
        return(input_fasta)

def run_blast(input,database,output,task,evalue,type = 'nucl',format = '10',container = None,threads = '1', debuglog = 'logs/debug_log.txt', repeatmasker = None):
    if type == 'nucl':
        blastype = 'blastn'
    elif type == 'prot':
        blastype = 'blastp'
    else:
        print('BLAST type not recognized, please specify "nucl" or "prot"')
        quit(1)
    with open(debuglog, 'a') as debug:
        if container is None:
            command = [blastype,'-db',database,'-query',input,'-out',output,'-outfmt',format,'-task',task,'-evalue',evalue,'-num_threads',threads]
            if repeatmasker is not None:
                command+=['-db_soft_mask','40']
            cdbg([command],debuglog)
        else:
            command = ['singularity','exec',container,blastype,'-db',database,'-query',input,'-out',output,'-outfmt',format,'-task',task,'-evalue',evalue]
            _ = debug.write(' '.join(command) + '\n')
            subprocess.call(command,stdout=debug, stderr=debug)
    print('Finished running BLAST!')
    print(f'Results file: {output}')
    print(f'Database location: {database}')
    print(f'Debug log and details: {debuglog}')

def type_check(inp):
    if inp in ['n','nucl','nucleotide']:
        outp = 'nucl'
    elif args.type in ['p','prot','protein']:
        outp = 'prot'
    else:
        print(f'Unrecognized BLAST search type "{args.type}"')
        quit(1)
    return(outp)

def task_check(inp,btype):
    nucl_options = ['dc-megablast','megablast','blastn','blastn-short']
    prot_options = ['blastp','blastp-short']
    if btype == 'default':
        if btype == 'nucl':
            outp = 'blastn'
        else:
            outp = 'blastp'
    else:
        if (inp in nucl_options and btype == 'nucl') or (inp in prot_options and bytpe == 'prot'):
            outp = inp
        else:
            print(f'BLAST task {inp} is not valid for search type {btype}')
    return(outp)


def mask_seq(input_fasta,bname,repeatlib,threads,container,debuglog):
    commandlist = []
    # get file name
    input_fasta_name = input_fasta.split('/')[-1]
    # use temporary directory for this process
    tempdir = f'temp/{bname}/'
    commandlist.append(['singularity','exec',container,'RepeatMasker','-xsmall','-dir',tempdir,'-lib',repeatlib,input_fasta,'-pa',threads])
    cdbg(commandlist,debuglog)
    # move masked fasta to new location
    #masked_fasta = temp_input_fasta + '.masked'
    masked_fasta = f'temp/{bname}/{input_fasta_name}.masked'
    output_fasta = input_fasta.replace('.fasta','_masked.fasta')
    if os.path.isfile(masked_fasta):
        # if the repeat masking worked, we can move the masked fasta back to db/ and generate the masking file
        commandlist2 = []
        # move the file back
        commandlist2.append(['mv',masked_fasta,output_fasta])
        # generate the masking file for makeblastdb
        mask_data = output_fasta.replace('.fasta','.asnb')
        commandlist2.append([
            'convert2blastmask','-in',output_fasta,'-masking_algorithm', 'repeat', '-masking_options', 
            'repeatmasker, default', '-outfmt', 'maskinfo_asn1_bin', '-out', mask_data])
        cdbg(commandlist2,debuglog)
        print(f'Masked fasta file created in {output_fasta}')
        print(f'Masking data created in {mask_data}')
    else:
        print( f'No repeats detected for {input_fasta}, the original file will be used for BLAST. Check the debug log to ensure that RepeatMasker detected the input file.')
        output_fasta = input_fasta
        mask_data = None
    return([output_fasta,mask_data])


def cdbg(clist,debuglog):
    with open(debuglog,'a') as debug:
        for command in clist:
            _ = debug.write(' '.join(command) + '\n')
            exitcode = subprocess.call(command,stdout=debug, stderr=debug)
            if exitcode != 0:
                print('Nonzero exit code when running the following command:')
                print(' '.join(command))
                print(f'Check debug log at {debuglog} for details')
                quit(1)

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
        default='NewSearch'
        )
    parser.add_argument(
        '--task','-k',type=str,choices=['default','dc-megablast','megablast','blastn','blastn-short','blastp','blastp-short'],
        help='''(Optional) Provide an alternate task for the BLAST search. The blastn-short/blastp-short options are for query sequences below 
        50 bases or 30 residues. The blastn/blastp options are for moderate-length queries. The megablast options are for very similar sequences, 
        with megablast for intraspecies searches and dc-megablast for cross-species. Both megablast options are only available for nucleotide BLAST.
        The default is blastn/blastp.''',
        default='default'
        )
    parser.add_argument(
        '--evalue','-e',type=str,help='''(Optional) Provide an evalue threshold. Hits above this threshold will not be displayed. 
        The default is 1e-3.''',
        default='1e-3'
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
        mask repetitive sequences in the query and subject. This is currently optimized for C. auris genomes. This will also enable 
        masking in an existing database.''',
        default=None
        )
    parser.add_argument(
        '--threads','-th',type=str,help='''(Optional) Provide the number of threads to use.''',
        default='1'
        )
    parser.add_argument(
        '--verbose', '-v', action='store_true', help='''(Optional) Enable verbose output.''',
        default=False
        )
    args = parser.parse_args()
    # handle any missing or incorrect arguments
    args.type = type_check(args.type)
    args.task = task_check(args.task,btype=args.type)
    if args.query is None or args.subject is None:
        print('Missing query or subject sequence')
        quit(1)
    # ensure query, subject, and annotation files exist, and convert them to absolute paths
    for f in [args.query,args.subject,args.annotation]:
        if f is not None:
            if not os.path.exists(f):
                print(f'Could not locate file or directory at {f}')
                quit(1)
    args.query,args.subject,args.annotation = [os.path.abspath(x) if x is not None else x for x in [args.query,args.subject,args.annotation]]
    # change working directory to location of script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # generate the files and directories needed for the BLAST search in this directory
    for d in ['logs/',f'results/{args.name}',f'db/{args.name}/',f'temp/{args.name}/']:
        if not os.path.isdir(d):
            subprocess.call(['mkdir','-p',d])
    # create debug log
    debuglog = f'logs/{args.name}_debuglog.txt'
    with open(debuglog,'w') as fh:
        _ = fh.write(f'debug log for {args.name}\n')
    # generate a single query file
    query_fasta = check_query(args.query,args.name,debuglog)
    if args.verbose:
        print(f'Query check complete, the query is currently {query_fasta}')
    # perform repeat masking on query, if needed
    repeatlib = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_11_11_funQCD_database/lib/repeat_libraries/fungi_b8441/b8441_fungi_repeatlib.fa'
    # repeat masking of the query does not appear to be compatible with command line BLAST
    # if args.repeatmasker is not None:
    #     query_fasta = mask_seq(
    #         input_fasta=query_fasta,bname=args.name,repeatlib=repeatlib,threads=args.threads,
    #         container=args.repeatmasker,debuglog=debuglog)
    #     if args.verbose:
    #         print(f'Query repeat masking complete. The query is currently {query_fasta}')
    # generate a BLAST database, if needed
    mask_data = None
    if not args.database:
        subject_fasta = check_subject(args.subject,args.name,debuglog)
        if args.verbose:
            print(f'Subject check complete. The subject is currently {subject_fasta}')
        if args.repeatmasker is not None:
            subject_fasta,mask_data = mask_seq(
                input_fasta=subject_fasta,bname=args.name,repeatlib=repeatlib,threads=args.threads,
                container=args.repeatmasker,debuglog=debuglog)
            if args.verbose:
                print(f'Subject repeat masking complete. The subject is currently {subject_fasta}')
        subject_fasta = make_blast_db(input_fasta = subject_fasta,output_name = args.name,type = args.type,container = args.container, debuglog = debuglog, mask_data = mask_data)
        if args.verbose:
            print(f'BLAST database generated. The subject is currently {subject_fasta}')
    else:
        subject_fasta = args.subject
        if args.verbose:
            print(f'Using an existing BLAST database. The subject is currently {subject_fasta}')
        if args.repeatmasker is not None:
            mask_data = 1
            if args.verbose:
                print('Enabling repeat masking in the existing database')
    search_output_name = f'results/{args.name}/{args.name}_{args.type}_{args.task}_{args.evalue}_blast_results.csv'
    run_blast(
        input=query_fasta, database=subject_fasta, output=search_output_name, task=args.task, type=args.type, evalue=args.evalue, 
        format=args.format, container=args.container, threads=args.threads, debuglog=debuglog, repeatmasker=mask_data)
    # addition needed here - attach annotation data to results

if __name__ == "__main__":
    main()

