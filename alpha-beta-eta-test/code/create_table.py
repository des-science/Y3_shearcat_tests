import fitsio
import argparse
import os
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Print table and save plot with it using samples of abe')

    parser.add_argument('--samplesabe', default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/sample_ab_Y3_03-31-20_JKmarco_2.fits')
    parser.add_argument('--nsig', default=1, type=int, 
                        help='How many sigmas for the marginalized confidence interval')
    parser.add_argument('--ndof', default=117, type=int, 
                        help='ndof just to put the label, 20*6-3')
    parser.add_argument('--burn_frac', default=0.05, type=float, 
                        help='Burn frac of samples, default 5%')
    parser.add_argument('--hartlap', default=1, type=float, 
                        help='hartlap factor for chi2')
    parser.add_argument('--outpath',
                        default='/data/git_repositories/Y3_shearcat_tests/alpha-beta-eta-test/measured_correlations/tomo_ab_margin_Y3_03-31-20_JKmarco/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--filename',
                        default='samples_Y3_03-31-20',
                        help='name of table without format (.tex and .png)')

    args = parser.parse_args()

    return args
 


def saveintex(models_combo, margin, overall, parlist, chisq_list, filename, ndof):
    print('Generating table.tex')
    eq, abe, ab, ae, be, a, b, e = models_combo
    if overall:
        parsbin1, parsbin2,  parsbin3,  parsbin4 =  parlist
        if abe:
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $\alpha$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$\\ \rule{0pt}{3ex} %s $\beta$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \rule{0pt}{3ex} %s $\eta$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline  %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n',  parsbin1[0], parsbin2[0], parsbin3[0], parsbin4[0],'\n',  parsbin1[1], parsbin2[1], parsbin3[1], parsbin4[1], '\n', parsbin1[2], parsbin2[2], parsbin3[2], parsbin4[2],'\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')
        if ab or ae or be:
            if ab: name1 = r'\alpha'; name2 = r'\beta'
            if ae: name1 = r'\alpha'; name2 = r'\eta'
            if be: name1 = r'\beta'; name2 = r'\eta'
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$\\ \rule{0pt}{3ex} %s $%s$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline  \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name1, parsbin1[0], parsbin2[0], parsbin3[0], parsbin4[0],'\n', name2,  parsbin1[1], parsbin2[1], parsbin3[1], parsbin4[1], '\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')
        if a or b or e:
            if a: name = r'\alpha'
            if b: name = r'\beta'
            if e: name = r'\eta'
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$\\  \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name, parsbin1[0], parsbin2[0], parsbin3[0], parsbin4[0],'\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')
        
        print(text[1:-1], file=open(filename, "w"))
        print(filename ,  'written!')
        
 

    if margin:
        parsbin1, parsbin2,  parsbin3,  parsbin4 =  parlist
        if abe:
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $\alpha$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$\\ \rule{0pt}{3ex} %s $\beta$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ \\ \rule{0pt}{3ex} %s $\eta$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ \\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n',  parsbin1[0][0],parsbin1[0][1],parsbin1[0][2],parsbin2[0][0],parsbin2[0][1],parsbin2[0][2],parsbin3[0][0],parsbin3[0][1],parsbin3[0][2],parsbin4[0][0],parsbin4[0][1],parsbin4[0][2],'\n',  parsbin1[1][0],parsbin1[1][1],parsbin1[1][2], parsbin2[1][0],parsbin2[1][1],parsbin2[1][2], parsbin3[1][0],parsbin3[1][1],parsbin3[1][2],parsbin4[1][0],parsbin4[1][1],parsbin4[1][2], '\n', parsbin1[2][0],parsbin1[2][1],parsbin1[2][2],parsbin2[2][0],parsbin2[2][1],parsbin2[2][2], parsbin3[2][0],parsbin3[2][1],parsbin3[2][2],parsbin4[2][0],parsbin4[2][1],parsbin4[2][2],'\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')
      
        if ab or ae or be:
            if ab: name1 = r'\alpha'; name2 = r'\beta'
            if ae: name1 = r'\alpha'; name2 = r'\eta'
            if be: name1 = r'\beta'; name2 = r'\eta'

            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$\\ \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ \\ \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name1,  parsbin1[0][0],parsbin1[0][1],parsbin1[0][2],parsbin2[0][0],parsbin2[0][1],parsbin2[0][2],parsbin3[0][0],parsbin3[0][1],parsbin3[0][2],parsbin4[0][0],parsbin4[0][1],parsbin4[0][2],'\n', name2,  parsbin1[1][0],parsbin1[1][1],parsbin1[1][2], parsbin2[1][0],parsbin2[1][1],parsbin2[1][2], parsbin3[1][0],parsbin3[1][1],parsbin3[1][2],parsbin4[1][0],parsbin4[1][1],parsbin4[1][2], '\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')

        if a or b or e:
            if a: name = r'\alpha'
            if b: name = r'\beta'
            if e: name = r'\eta'

            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|c|c|c|} \hline %s & \textrm{Bin}1 & \textrm{Bin}2 & \textrm{Bin}3 & \textrm{Bin}4 \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$ & $%.3f_{%.3f}^{+%.3f}$\\ \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name,  parsbin1[0][0],parsbin1[0][1],parsbin1[0][2],parsbin2[0][0],parsbin2[0][1],parsbin2[0][2],parsbin3[0][0],parsbin3[0][1],parsbin3[0][2],parsbin4[0][0],parsbin4[0][1],parsbin4[0][2],'\n', ndof*chisq_list[0],  ndof*chisq_list[1], ndof*chisq_list[2], ndof*chisq_list[3], '\n', ndof, chisq_list[0],  chisq_list[1], chisq_list[2], chisq_list[3],'\n', '\n')

        print(text[1:-1], file=open(filename, "w"))
        print(filename ,  'written!')

def saveintex_nontomo(models_combo, margin, overall, parlist, chisq_list, filename, ndof):
    print('Generating table.tex')
    eq, abe, ab, ae, be, a, b, e = models_combo
    if overall:
        if abe:
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic} \\ \hline \rule{0pt}{3ex} %s $\alpha$ & $%.3f$ \\ \rule{0pt}{3ex} %s $\beta$ & $%.3f$ \\ \rule{0pt}{3ex} %s $\eta$ & $%.3f$ \\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$\\ \hline  %s $\chi^{2}_{\nu=%d}$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n',  parlist[0],  '\n',  parlist[1],  '\n', parlist[2], '\n', ndof*chisq_list, '\n', ndof, chisq_list,'\n', '\n')
        if ab or ae or be:
            if ab: name1 = r'\alpha'; name2 = r'\beta'
            if ae: name1 = r'\alpha'; name2 = r'\eta'
            if be: name1 = r'\beta'; name2 = r'\eta'
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic}  \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f$ \\ \rule{0pt}{3ex} %s $%s$ & $%.3f$\\ \hline  \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name1, parlist[0],'\n', name2,  parlist[1],  '\n', ndof*chisq_list,  '\n', ndof, chisq_list,'\n', '\n')
        if a or b or e:
            if a: name = r'\alpha'
            if b: name = r'\beta'
            if e: name = r'\eta'
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic}  \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f$\\  \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name, parlist[0],'\n', ndof*chisq_list,  '\n', ndof, chisq_list, '\n', '\n')
        
        print(text[1:-1], file=open(filename, "w"))
        print(filename ,  'written!')
        
 

    if margin:
        if abe:
            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic}  \\ \hline \rule{0pt}{3ex} %s $\alpha$ & $%.3f_{%.3f}^{+%.3f}$\\ \rule{0pt}{3ex} %s $\beta$ & $%.3f_{%.3f}^{+%.3f}$ \\ \rule{0pt}{3ex} %s $\eta$ & $%.3f_{%.3f}^{+%.3f}$\\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$\\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n',  parlist[0][0],parlist[0][1],parlist[0][2],'\n',  parlist[1][0],parlist[1][1],parlist[1][2], '\n', parlist[2][0],parlist[2][1],parlist[2][2], '\n', ndof*chisq_list,  '\n', ndof, chisq_list,'\n', '\n')
      
        if ab or ae or be:
            if ab: name1 = r'\alpha'; name2 = r'\beta'
            if ae: name1 = r'\alpha'; name2 = r'\eta'
            if be: name1 = r'\beta'; name2 = r'\eta'

            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic}  \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$ \\ \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$\\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$ \\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name1,  parlist[0][0],parlist[0][1],parlist[0][2], '\n', name2,  parlist[1][0],parlist[1][1],parlist[1][2], '\n', ndof*chisq_list , '\n', ndof, chisq_list, '\n', '\n')

        if a or b or e:
            if a: name = r'\alpha'
            if b: name = r'\beta'
            if e: name = r'\eta'

            text =  r"$\begin{center} %s \centering %s \begin{tabular}{ |c|c|} \hline %s & \textrm{Non-tomographic}   \\ \hline \rule{0pt}{3ex} %s $%s$ & $%.3f_{%.3f}^{+%.3f}$\\ \hline \rule{0pt}{3ex} %s $\chi^{2}$ & $%.3f$ \\ \hline %s $\chi^{2}_{\nu=%d}$ & $%.3f$\\ \hline %s \end{tabular} %s \end{center}$"%('\n','\n','\n', '\n', name,  parlist[0][0],parlist[0][1],parlist[0][2], '\n',  ndof*chisq_list, '\n', ndof, chisq_list, '\n', '\n')

        print(text[1:-1], file=open(filename, "w"))
        print(filename ,  'written!')

def pdflatex_standalone(bodytexfile,  outpath, filename):
    import subprocess
    import datetime
    import shutil
    ##create auxdirectory
    starttime = datetime.datetime.now()
    prefix = starttime.strftime("%Y%m%dT%H%M%S")
    auxpath = os.path.join(outpath, 'latex_tempfiles_%s'%(prefix))
    subprocess.call(['mkdir', auxpath])
    
    #create tex file compilable with pdflatex
    main_tex = "%s_aux.tex"%(filename)
    template = r'''\documentclass[preview]{{standalone}}
\usepackage{{booktabs}}
\begin{{document}}
\centering
{}
\end{{document}}
'''.format(bodytexfile)
    print(template[:], file=open(os.path.join(auxpath, main_tex), "w"))

    ##saving all latex file in auxpath, and convert pdf to png
    subprocess.call(['pdflatex', '-output-directory',auxpath,'-jobname', filename, main_tex])
    print("PDFLATEX FINISHED")
    pdffile = os.path.join(auxpath, "%s.pdf"%(filename))
    outname = os.path.join(outpath, "%s.png"%(filename))
    subprocess.call(['convert', '-density','300', pdffile, '-strip', outname])
    print("CONVERT FINISHED")
    #remove temporary folder
    shutil.rmtree(auxpath)

def main():
    from src.chi2 import chi2nu,  ndof
    from src.readfits import  read_rhos, read_taus
    args = parse_args()

    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    abesamps = fitsio.read(args.samplesabe)
    df = pd.DataFrame(abesamps)
    print(df)
    #print(df.columns)
    #Chi2 info should be save

    i = 0
    while True:
        i += 1
        if 'chi2nu%i'%(i) not in abesamps.dtype.names: break
    nbins = i - 1
    
    af,bf,ef, xf = [name in abesamps.dtype.names for name in ['a1', 'b1', 'e1', 'chi2nu1'] ]
    ext = '%s%s%s'%('a'*af,'b'*bf, 'e'*ef)

    tablename = '%s%s'%(args.filename,'.tex')
    plotname = '%s%s'%(args.filename,'.png')

    if len(abesamps.dtype.names) > 7: bnames = ['Bin %i'%(i) for i in range(1, nbins + 1)]
    else: bnames = ['Non-tomographic']
    cols = ["Parameter"] +  bnames
    print(cols)
   
    if args.nsig == 1: percents = [16, 50, 84]
    elif args.nsig ==2: percents = [2.3, 50, 97.7]
    else: print("Not defined percentiles")

    res = '%.3f'
    alphas = af*[r"$\alpha$"]
    betas = bf*[r"$\beta$"]
    etas = ef*[r"$\eta$"]
    chis = xf*[r"$\chi^{2}_{\nu=%i}$"%(args.ndof)]
    if af:
        for col in ['a%i'%(i) for i in range(1, nbins + 1)]:
            l, m, r = np.percentile(abesamps[col], percents)
            alphas += [r'$%s_{-%s}^{+%s}$'%(res,res,res)%(m,m - l,r - m)]
    if bf:
        for col in ['b%i'%(i) for i in range(1, nbins + 1)]:
            l, m, r = np.percentile(abesamps[col], percents)
            betas += [r'$%s_{-%s}^{+%s}$'%(res,res,res)%(m,m - l,r - m)]
    if ef:
        for col in ['e%i'%(i) for i in range(1, nbins + 1)]:
            l, m, r = np.percentile(abesamps[col], percents)
            etas += [r'$%s_{-%s}^{+%s}$'%(res,res,res)%(m,m - l,r - m)]
    if xf:    
        for col in ['chi2nu%i'%(i) for i in range(1, nbins + 1)]:
            l, m, r = np.percentile(args.hartlap*abesamps[col], percents)
            chis += [r'$%s_{-%s}^{+%s}$'%(res,res,res)%(m,m - l,r - m)]

    #print("Alphas :\n",  alphas)
    #print("Betas :\n",  betas)
    #print("Etas :\n",  etas)
    #print("chis :\n",  chis)

    
    
    df = pd.DataFrame([x for x in [alphas, betas, etas, chis] if x], columns=cols)
    #print(df)
    tabletex = df.to_latex(index=False, escape=False)
    print(tabletex) 

    print(tabletex[:], file=open(os.path.join(outpath, "%s.tex"%(args.filename)), "w"))

    pdflatex_standalone(tabletex,  outpath, args.filename)
    
     
    
    
 

    
if __name__ == "__main__":
    main()
