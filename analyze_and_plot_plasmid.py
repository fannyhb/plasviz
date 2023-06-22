from read_fasta import read_fasta

import subprocess as sp
import shlex
import argparse
from os import path
from pycirclize import Circos
from Bio.SeqFeature import SeqFeature, FeatureLocation
from matplotlib.patches import Patch


isfinder_db = '/storage/shared/databases/isfinder/ISFinder.faa.dmnd'
nr_db = '/storage/shared/databases/nr_DB/ncbi_protein_bacterial_c95.faa.dmnd'
mobileOG_db = '/storage/fannyb/databases/mobileOG/man_cur/mobileOG_plasmid_build.dmnd'
resfinder_db = '/storage/fannyb/databases/resfinder/2022_11/resfinder_inti1/pure_resfinder_plasmid_build.dmnd'
scov = 60
base_cutoff = 60
threads = 30
shift = True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',required=True)
    parser.add_argument('-o','--outdir',default='./')
    

    options = parser.parse_args()

    outdir = path.abspath(options.outdir)
    infile = options.infile
    basename = '.'.join(options.infile.split('/')[-1].split('.')[:-1])
    orffile = f'{outdir}/{basename}_orfs.fasta'
    prodigal_listfile = f'{outdir}/{basename}.prodigal'

    orf_info = {}
    run_prodigal(infile,orffile,prodigal_listfile)
    orf_info, seqlen = parse_prodigal(orffile, orf_info, prodigal_listfile)

    if shift:
        new_fasta = f'{outdir}/{basename}_shifted.fasta'
        shift_len(orf_info, infile, new_fasta)
        orf_info = {}
        run_prodigal(new_fasta,orffile,prodigal_listfile)
        orf_info, seqlen = parse_prodigal(orffile, orf_info, prodigal_listfile)

    dbs = {'isfinder':isfinder_db,'nr':nr_db,'mobileOG':mobileOG_db,'resfinder':resfinder_db}

    for db, db_path in dbs.items():
        outfile = f'{outdir}/{basename}.{db}'
        if not path.isfile(outfile):
            run_diamond(orffile,db_path,base_cutoff,scov,threads,outfile)
        orf_info = parse_diamond(outfile, 70,db, orf_info)
    for orf, info in orf_info.items():
        print(orf, info)
    
    outfig = f'{outdir}/{basename}.png'
    plot_the_plasmid(orf_info, seqlen, basename,outfig)

def plot_the_plasmid(orf_info, seqlen, basename, outfig):
    db_colors = {'resfinder':"firebrick",'isfinder': "goldenrod",'mobileOG':"olivedrab",'nr':"cornflowerblue"}
    text_inside = False
    sectors = {'my_plasmid': seqlen}
    circos = Circos(sectors)
    circos.text(f'{basename}: {seqlen}bp',size=12,r=20)
    
    sector = circos.sectors[0]
    cds_track = sector.add_track((90,100))
    cds_track.axis(ec="none")

    plotted = []
    text_pos = []
    labels = []
    for db, arrow_col in db_colors.items():
        cds_track, plotted, pos, names = plot_resfinder(orf_info, cds_track, db, arrow_col,plotted,text_inside)
        text_pos.extend(pos)
        labels.extend(names)
    
# Plot xticks & intervals on inner position
    cds_track.xticks_by_interval(
        interval=1000,
        outer=False,
        show_bottom_line=True,
        label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
        )
    
    if not text_inside:
        # Plot CDS product labels on outer position
        cds_track.xticks(
            text_pos,
            labels,
            label_orientation="horizontal",
            show_bottom_line=False,
            label_size=6,
            line_kws=dict(ec="grey"),
    )
    
    fig = circos.plotfig()
    # add legends
    handles = [
            Patch(color=db_colors['resfinder'],label='ARG'),
            Patch(color=db_colors['isfinder'],label='IS'),
            Patch(color=db_colors['mobileOG'],label='MGE'),
            Patch(color=db_colors['nr'],label='Misc protein')
            ]
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=11)

    fig.savefig(outfig)

def shift_len(orf_info, infile, outfile):
    fout = open(outfile,'w')
    starts = []
    for orf, info in orf_info.items():
        if info['complete']:
            starts.append(info['start'])
    start = min(starts)          
#    start = max(starts)          
    for header, seq in read_fasta(infile,False):
        new_start = start-1
        start_seq = seq[0:new_start]
        end_seq = seq[new_start:] + start_seq
        fout.write(f'>{header}\n{end_seq}\n')

def plot_resfinder(orf_info, cds_track,db,arrow_col,plotted,inside):
    cds_list = []
    hit = False
    pos_list, name_list = [], []
    for orf, info in orf_info.items():
        if db in info:
            if not orf in plotted:
                cds_list.append(SeqFeature(FeatureLocation(info['start'],info['end'],strand=info['strand']),type="CDS"))
                plotted.append(orf)
                hit = True
                text_pos = info['start'] + round((info['end']-info['start'])/2)
                if inside:
                    cds_track.text(info[db]['name'],text_pos, size=8,color="black")
                else:
                    pos_list.append(text_pos)
                    name_list.append(info[db]['name'])
    if hit:
        cds_track.genomic_features(cds_list,
                plotstyle="arrow",
                r_lim=(90,100),
                fc=arrow_col,
                )
    return cds_track, plotted, pos_list, name_list

def parse_diamond(infile,cutoff,db,orf_info):
    orf_gene_info = {}
    with open(infile) as f:
        for line in f:
            populate = False
            line = line.strip().split('\t')
            orf = line[0]
            if db == 'resfinder':
                gene = line[1]
                accession = gene.split('_')[-1]
                name = gene.split('_')[0]
            elif db == 'mobileOG':
                gene = line[1]
                name = gene.split('|')[1]
                accession = gene.split('|')[2]
            elif db == 'nr':
                gene = line[1]
                accession = line[1]
                name = '_'.join(line[-1].split('[')[0].split('(')[0].strip().split())
                name = '_'.join(name.split('_')[1:])
                if 'hypothetical_protein' in name:
                    name = 'hypothetical_protein'
                    name = 'hyp. prot.'
                elif name.startswith('type'):
                    name = '_'.join(name.split('_')[0:3])
                elif name.startswith('plasmid'):
                    name = '_'.join(name.split('_')[0:3])
                else:
                    name = name.split('_')[0]
            elif db == 'isfinder':
                gene = line[1]
                accession = line[1]
                name = line[1]
            else:
                gene = line[1]
                accession = gene.split('_')[-1]
                name = gene.split('_')[0]
            pid = line[2]
            aln = line[3]
            score = float(line[11])
            if not orf in orf_gene_info:
                gene_info = {}
                populate = True
            else:
                gene_info = orf_gene_info[orf]
                if score > float(gene_info['score']):
                    populate = True
            if populate:
                gene_info['pid'] = float(pid)
                gene_info['aln'] = int(aln)
                gene_info['gene'] = gene
                gene_info['name'] = name
                gene_info['accession'] = accession
                gene_info['score'] = score
                orf_gene_info[orf] = gene_info
                print(f'pop orf: {orf}')
    # Move my new dictionary to the master dictionary
    for orf, info in orf_info.items():
        if orf in orf_gene_info:
            orf_info[orf][db] = orf_gene_info[orf]
    return orf_info

def parse_prodigal(infile,orf_info, prod_list):
    for header, seq in read_fasta(infile,False):
        line = header.split()
        orf = line[0]
        tmp_dict = {}
        tmp_dict['start'] = int(line[2])
        tmp_dict['end'] = int(line[4])
        tmp_dict['strand'] = int(line[6])
        tmp_dict['length'] = tmp_dict['end']-tmp_dict['start']+1
        tmp = line[8].split(';')[1].split('=')[-1]
        if tmp == '00':
            full=True
        else:
            full=False
        tmp_dict['complete'] = full
        orf_info[orf] = tmp_dict
    with open(prod_list) as f:
        for line in f:
            if line.startswith('DEFINITION'):
                line = line.split()[1].strip().split(';')
                seqlen = int(line[1].split('=')[-1])                
    return orf_info, seqlen

def run_diamond(infile,database,identity,scov,threads,outfile):
    msg = f'diamond blastx --db {database} -q {infile} -o {outfile}'\
            f' --id {identity} --subject-cover {scov} --threads {threads}'\
            f' --outfmt 6 qseqid sseqid pident length mismatch gapopen'\
            f' qstart qend sstart send evalue bitscore qlen slen qframe'\
            f' qtitle stitle'
    sp.run(shlex.split(msg))

def run_prodigal(infile,nuc_outfile,list_outfile):
    msg = f'prodigal -i {infile} -o {list_outfile}'\
            f' -p meta -d {nuc_outfile}'
    sp.run(shlex.split(msg))

if __name__=='__main__':
    main()
