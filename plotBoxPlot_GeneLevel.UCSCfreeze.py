
import scipy as sp
import pdb
import libs.expression

import os
import sys

from optparse import OptionParser, OptionGroup


def parse_options(argv):
    '''
    option parser
    '''
    ### required input
    parser = OptionParser()
    required  = OptionGroup(parser, 'Input')
    required.add_option('-s', '--star', dest = 'fn_star', metavar = 'FILE', help = 'STAR count file') 
    required.add_option('-t', '--tophat', dest = 'fn_tophat', metavar = 'FILE', help = 'Tophat count file') 


    ### optional arguments
    optional = OptionGroup(parser, 'optional')
    
    ## outdir
    outdir_hstr = 'Output direcory for all files'
    outdir = os.path.realpath(__file__).rsplit('/',1)[0]
    optional.add_option('-o', '--out_dir', dest = 'out_dir', help = outdir_hstr, default = outdir)

    ## genelength buffe
    fnl_hstr = 'Location of file with gene length'
    fnlfile =  os.path.join(os.path.realpath(__file__).rsplit('/',1)[0], 'data','geneLength.tsv')
    optional.add_option('-g', '--genelength', dest = 'fn_length',metavar = 'FILE', help = fnl_hstr, default = fnlfile)

    ## annotation file
    fna_hstr = 'Annotation file [gtf]'
    fnafile = os.path.join(os.path.realpath(__file__).rsplit('/',1)[0], 'data','gencode.v19.annotation.hs37d5_chr.gtf')
    optional.add_option('-a', '--anno', dest = 'fn_anno', metavar = 'FILE', help = fna_hstr, default = fnafile)

    
    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    
    return options

def getLength(fn_anno):
    ''' 
    input gtf annotation and output all exonic positions
    '''

    anno  = sp.loadtxt(fn_anno, delimiter = '\t', dtype = 'string', usecols=[0,2,3,4,8])
    chrm = anno[:,0]
    anno = anno[:,1:]
    anno  = anno[anno[:,0] == 'exon',:] ### filter down to exons
    pos   = anno[:,1:3].astype('int')
    
    gid   = [x.split(';')[0] for x in anno[:,3]] ### clean gene id's
    gid   = sp.array([x.split(" ")[1].strip('\"') for x in gid])
    
    ugid  = sp.unique(gid)
    gnL   = sp.zeros(ugid.shape[0])
    for i,mgid in enumerate(ugid):
        if (i % 100) == 0:
            sys.stdout.write("%i out of %i \n" % (i, ugid.shape[0]))
        mpos   = sp.hstack([range(x[0],x[1]+1) for x in pos[mgid == gid,:]]) ### positions for this frame
        gnL[i] = sp.unique(mpos).shape[0]
    return ugid, gnL

def getSizeFactor(fn_anno, data, gid, mode = 'sum', withXYMT = True, filterbyPC = True):
    '''
    input annotation, counts and gene ids
    output sum of protein coding gene levels excluding sex chromosomes and mitochondria genes
    '''
    anno  = sp.loadtxt(fn_anno, delimiter = '\t', dtype = 'string', usecols=[0,2,8])
    anno  = anno[anno[:,1] == 'gene', :]
    if not withXYMT: ### filter xymt
        anno  = anno[anno[:,0] != 'MT',:]
        anno  = anno[anno[:,0] != 'Y',:]
        anno  = anno[anno[:,0] != 'X',:]

    agid   = [x.split(';')[0] for x in anno[:,2]] ### clean gene id's
    agid   = sp.array([x.split(" ")[1].strip('\"') for x in agid])

    if filterbyPC: ### filter protein coding
        gtpe  = [x.split(';')[2] for x in anno[:,2]]
        gtpe  = sp.array([x.split('\"')[1].split('\"')[0] for x in gtpe])
        iPC   = sp.where(gtpe == 'protein_coding')[0]
        agid  = agid[iPC]

    iGn = sp.in1d(gid, agid)
    libsize = sp.sum(data[iGn,:], axis = 0) 
    if mode == 'uq':
         libsize = sp.array([sp.percentile(x[x!=0] ,75) for x in data[iGn,:].T])  * iGn.sum() 

    return libsize



def main():
    ### get options
    options = parse_options(sys.argv)

    ### get gene length and store locally (Should probably clean itself up)
    if os.path.exists(options.fn_length):
        tmp = sp.loadtxt(options.fn_length, delimiter = '\t', dtype = 'string')
        gnLength     = tmp[:,1].astype('float')
        gnLength_gid = tmp[:,0]
    else:
        gnLength_gid, gnLength = getLength(options.fn_anno)
        sp.savetxt(options.fn_length, sp.vstack((gnLength_gid, gnLength)).T, delimiter = '\t', fmt = '%s')


    ### get star data
    data = sp.loadtxt(options.fn_star, delimiter = '\t', dtype = 'string')
    gtid = data[0,1:]
    gid  = data[1:,0]
    data = data[1:,1:].astype('float')

    ### remove htseq summary stats if available
    iOK  = ~sp.array([x.startswith('__') for x in gid])
    gid  = gid[iOK]
    data = data[iOK,:]

    ### sort by gene id
    sidx = sp.argsort(gid)
    gid  = gid[sidx]   
    data = data[sidx,:]

    ### sort by gtid
    sidx = sp.argsort(gtid)
    gtid = gtid[sidx]
    data = data[:,sidx]

    ### keep copy for later and sort gene length by gene id so that id's match 
    ### in case length is ordered differently
    data_star    = data.copy()
    gtid_star    = gtid.copy()
    sidx         = sp.argsort(gnLength_gid)
    gnLength_gid = gnLength_gid[sidx]
    gnLength     = gnLength[sidx]

    ### get total count FPKM
    libSize  = getSizeFactor(options.fn_anno, data, gid, withXYMT = True, filterbyPC = False)     

    ### store size factor total count star
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_star_gene_fpkm'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'star_fpkm.tsv.gz'))

    ### get uq FPKM
    libSize = getSizeFactor(options.fn_anno, data, gid, mode = 'uq', withXYMT = False, filterbyPC = True)      ###'uq'
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_star_gene_fpkm_uq'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'star_fpkm_uq.tsv.gz'))

    data = sp.loadtxt(options.fn_tophat, delimiter = '\t', dtype = 'string')
    gtid = data[0,1:]
    gid  = data[1:,0]
    gid  = sp.array([x.strip('\"') for x in gid])

    data = data[1:,1:].astype('float')
    iOK  = ~sp.array([x.startswith('__') for x in gid])
    gid  = gid[iOK]
    data = data[iOK,:]
    sidx = sp.argsort(gid)
    gid  = gid[sidx]
    data = data[sidx,:]

    sidx = sp.argsort(gtid)
    gtid = gtid[sidx]
    data = data[:,sidx]


    if gtid_star.shape[0] != gtid.shape[0]:
        midx = sp.in1d(gtid, gtid_star)
        gtid = gtid[midx]
        data = data[:,midx]


    ### get size factor TOTAL COUNT based on PROTEin CODinG only
    libSize = getSizeFactor(options.fn_anno, data, gid, withXYMT = True, filterbyPC = False) 
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_th_gene_fpkm'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'th_fpkm.tsv.gz'))


    libSize = getSizeFactor(options.fn_anno, data, gid, mode = 'uq', withXYMT = False, filterbyPC = True)     
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_th_gene_fpkm_uq'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'th_fpkm_uq.tsv.gz'))
    

    ### find overlap
    instar = sp.in1d(gtid, gtid_star)
    intophat = sp.in1d(gtid_star, gtid)
    assert instar.sum() == intophat.sum(), 'Samples do not match'
    assert instar.sum() == gtid.shape[0], 'Partial mismatch of ids'
    assert intophat.sum() == gtid_star.shape[0], 'Partial mismatch of ids'


    ### subset to intersection
    data = data[:,instar]
    data_star = data_star[:,intophat]
    
    ### fix labels
    gtid = gtid[instar]
    gtid_star = gtid_star[intophat]
    ### done matching
    data = data + data_star
    libSize = getSizeFactor(options.fn_anno, data, gid, withXYMT = True, filterbyPC = False) 
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_joint_gene_fpkm'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'joint_fpkm.tsv.gz'))

    libSize  = getSizeFactor(options.fn_anno, data, gid, mode = 'uq', withXYMT = False, filterbyPC = True)     
    myExp = libs.expression.ExpressionData(Y = data, GID = gid, GTID = gtid)
    myExp.libsizenorm(trafo='fpkm', totalcounts = libSize, gelen =gnLength)
    myExp.libsizeplots(fn_base = os.path.join(options.out_dir,'postnorm_joint_gene_fpkm_uq'), fmt = 'pdf', figsize = [100,5])
    myExp.store(os.path.join(options.out_dir, 'joint_fpkm_uq.tsv.gz'))



if __name__ == "__main__":
    main()
