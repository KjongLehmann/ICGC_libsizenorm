import scipy as sp
import os
import pdb

class ExpressionData():
    """
    A class for expression level data
    """


    def __init__(self, Y = None, GID = None, GTID = None):
        """
        Y = genes x samples
        """

        ## TODO: include a assert for dimensiosn such that Y is correct order
        self.Y    = Y
        self.GID  = GID
        self.GTID = GTID
    
    def store(self,fn):
        sp.savetxt(fn, sp.hstack((self.GID[:,sp.newaxis], self.Y)) , header = "\t".join(sp.hstack(("feature", self.GTID))) , comments = '',  delimiter = '\t', fmt = '%s')

    ### lib size normalization functionalities
    def libsizenorm(self, trafo = 'uq', totalcounts = None, gelen = None, return_libsize=False):
        if trafo == 'fpkm':
            assert totalcounts is not None, 'Error, Fpkm requires totalcounts'
            assert gelen is not None, 'Error, FPKM requires gene/exon length'
            libSize = totalcounts
            self.Y = (self.Y * 1E9) / (totalcounts)
            self.Y = (self.Y.T / gelen).T            
        if return_libsize:
            return (self.Y, totalcounts)

    def libsizeplots(self, fn_base = 'libsize', trafo = 'log',fmt = 'pdf', figsize = None, ylim = None, sort = True, ylab='Log2 Normalized Expression', xlab='Samples', colorlabels = None):
        ### import plotting libs
        import matplotlib
        matplotlib.use('AGG')
        import matplotlib.pyplot as plt

        if trafo == 'log':
            Y_trafo = sp.log2(self.Y + 1)
        else:
            Y_trafo = self.Y

        if figsize is None:
            fig = plt.figure()
        else:
            fig = plt.figure(figsize = figsize)
        
        if sort:
            iqr = sp.array([sp.median(x) for x in Y_trafo.T])
            sidx = sp.argsort(iqr)
            Y_trafo = Y_trafo[:,sidx]
        ax  = fig.add_subplot(111)
        if colorlabels is None:
            ax.boxplot([x for x in Y_trafo.T],sym = '', showfliers = False)
        else:
            colorlabeluq = sp.unique(colorlabels)
            ax.set_xlim(0,colorlabels.shape[0]+1)
            mymap = plt.get_cmap('ocean')
            for i,label in enumerate(colorlabeluq):
                midx =label == colorlabels       
                boxprops = dict(color = mymap(i/float(colorlabeluq.shape[0])))
                ax.boxplot([x for x in Y_trafo.T[midx,:]],sym = '', showfliers = False,positions=sp.where(midx)[0]+1, boxprops = boxprops, showmeans = True)

        if not ylim is None:
            ax.set_ylim(ylim)
        ax.set_yticklabels(ax.get_yticks().astype('|S10'))
        ax.set_ylabel(ylab)
        ax.set_xticklabels([])
        fn_base = fn_base + '_boxplot.' + fmt
        plt.savefig(fn_base, format = fmt)


